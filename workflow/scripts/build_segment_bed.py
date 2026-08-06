"""Build the segment BED shared by all bulk assays.

Runpeng Luo (2026-07-22)

Starts from the accessible-regions BED (region.bed; its 4th column is the
chromosome-arm ``region_id``), subtracts the blacklist, then splits each arm at the
union of all datasets' BEDPE breakpoints. Two ids are emitted per segment:

  region_id : the arm id carried from the input BED (RDR/QC grouping bound).
  seg_id    : ``{region_id}#{k}``, the breakpoint chunk (the hard bin boundary);
              k increments at each breakpoint cut inside the arm. With no BEDPE,
              every arm has one segment, so ``seg_id`` == one-per-arm (current
              behavior). BEDPE and BED share 0-based coordinates.

Dependencies:
  pandas, numpy; io_utils.read_BED, io_utils.read_BEDPE, utils.setup_logging.

Inputs
  region_bed: accessible-regions BED (BED4; 4th column = arm region_id).
  blacklist_bed: optional BED of regions to subtract.
  bedpe: zero or more BEDPE files of SV breakpoints.
Outputs:
  segment_bed: BED5 (#CHR, START, END, region_id, seg_id).
"""

import logging
import os

import numpy as np
import pandas as pd

from io_utils import read_BED, read_BEDPE
from utils import setup_logging

snakemake_handle = snakemake  # noqa: F821
setup_logging(snakemake_handle.log[0])

region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = snakemake_handle.input.get("blacklist_bed", None)
bedpe_files = list(snakemake_handle.input.get("bedpe", []))
out_bed = snakemake_handle.output["segment_bed"]


def subtract_blacklist(regions: pd.DataFrame, bl: pd.DataFrame) -> pd.DataFrame:
    """Interval-subtract blacklist per chromosome, preserving region_id."""
    bl_by_chr = {c: g[["START", "END"]].to_numpy() for c, g in bl.groupby("#CHR")}
    rows = []
    for chrom, start, end, rid in zip(
        regions["#CHR"].astype(str),
        regions["START"].astype(int),
        regions["END"].astype(int),
        regions["region_id"],
    ):
        segments = [(int(start), int(end))]
        for bs, be in bl_by_chr.get(chrom, ()):
            nxt = []
            for s, e in segments:
                if be <= s or bs >= e:
                    nxt.append((s, e))
                    continue
                if bs > s:
                    nxt.append((s, int(bs)))
                if be < e:
                    nxt.append((int(be), e))
            segments = nxt
        for s, e in segments:
            if e > s:
                rows.append([chrom, s, e, rid])
    return pd.DataFrame(rows, columns=["#CHR", "START", "END", "region_id"])


# arm regions, chr-normalized (region_id = 4th column of the accessible-regions BED)
regions = read_BED(region_bed)[["#CHR", "START", "END", "region_id"]].copy()
n_raw = len(regions)

if blacklist_bed:
    bl = read_BED(blacklist_bed, extra_columns=())[["#CHR", "START", "END"]]
    regions = subtract_blacklist(regions, bl)
    logging.info(
        f"blacklist subtracted: {n_raw} -> {len(regions)} region pieces "
        f"({os.path.basename(blacklist_bed)})"
    )

# union of per-chromosome breakpoint cut positions (each BEDPE end contributes a cut)
cuts = {}
for path in bedpe_files:
    bp = read_BEDPE(path)
    n_bp = 0
    for chrom_col, start_col in (("#CHR1", "START1"), ("#CHR2", "START2")):
        if chrom_col in bp.columns:
            for chrom, pos in zip(bp[chrom_col].astype(str), bp[start_col].astype(int)):
                cuts.setdefault(chrom, set()).add(int(pos))
                n_bp += 1
    logging.info(
        f"{os.path.basename(path)}: {len(bp)} junctions, {n_bp} breakpoint ends"
    )

cuts = {c: np.array(sorted(v), dtype=np.int64) for c, v in cuts.items()}
logging.info(
    f"union: {sum(len(v) for v in cuts.values())} breakpoints across {len(cuts)} chromosomes"
)

# split each region piece at cuts inside it; seg index counts cuts at or before the
# piece start, so all pieces of one arm in the same inter-cut range share a seg_id
rows = []
for chrom, start, end, rid in zip(
    regions["#CHR"].astype(str),
    regions["START"].astype(int),
    regions["END"].astype(int),
    regions["region_id"],
):
    chrom_cuts = cuts.get(chrom, np.empty(0, dtype=np.int64))
    inside = [int(c) for c in chrom_cuts if start < c < end]
    bounds = [int(start), *inside, int(end)]
    for a, b in zip(bounds[:-1], bounds[1:]):
        seg_k = int(np.searchsorted(chrom_cuts, a, side="right"))
        rows.append([chrom, a, b, rid, f"{rid}#{seg_k}"])

out = pd.DataFrame(rows, columns=["#CHR", "START", "END", "region_id", "seg_id"])
out = out.sort_values(["#CHR", "START"]).reset_index(drop=True)

seg_len = out["END"] - out["START"]
logging.info(
    f"segment bed: {len(regions)} region pieces -> {len(out)} segments, "
    f"{out['seg_id'].nunique()} seg_ids, {out['region_id'].nunique()} region_ids; "
    f"segment length bp: min={int(seg_len.min())}, median={int(seg_len.median())}, "
    f"mean={seg_len.mean():.0f}, max={int(seg_len.max())}"
)

out.to_csv(out_bed, sep="\t", header=False, index=False)
