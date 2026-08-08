"""Build the segment BED shared by every workflow mode.

Runpeng Luo
Last update: 2026-08-07

Takes the configured segmentation (``segment_bed``; its 4th column is the ``seg_id``),
stamps the chromosome arm each segment sits in, and subtracts the blacklist. Two ids are
emitted per segment:

  region_id : the arm id from region_bed (RDR/QC grouping bound).
  seg_id    : carried from segment_bed; the hard bin boundary, binning never merges a bb
              across it. With segment_bed == region_bed there is one segment per arm and
              seg_id == region_id.

Blacklist subtraction may split one segment into several rows; they all keep its seg_id,
so a blacklist hole never becomes a bin boundary.

Dependencies:
  pandas; io_utils.read_BED, range_utils.assign_range_to_range, utils.setup_logging.

Inputs
  segments: the configured segmentation BED (BED4; 4th column = seg_id).
  region_bed: accessible-regions BED (BED4; 4th column = arm region_id).
  blacklist_bed: optional BED of regions to subtract.
Outputs:
  segment_bed: BED5 (#CHR, START, END, region_id, seg_id).
"""

import logging
import os

import pandas as pd

from io_utils import read_BED
from range_utils import assign_range_to_range
from utils import setup_logging

snakemake_handle = snakemake  # noqa: F821
setup_logging(snakemake_handle.log[0])

segments_bed = snakemake_handle.input["segments"]
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = snakemake_handle.input.get("blacklist_bed", None)
out_bed = snakemake_handle.output["segment_bed"]

ID_COLS = ["region_id", "seg_id"]


def subtract_blacklist(regions: pd.DataFrame, bl: pd.DataFrame) -> pd.DataFrame:
    """Interval-subtract blacklist per chromosome, preserving ID_COLS."""
    bl_by_chr = {c: g[["START", "END"]].to_numpy() for c, g in bl.groupby("#CHR")}
    rows = []
    for chrom, start, end, *ids in regions[
        ["#CHR", "START", "END", *ID_COLS]
    ].itertuples(index=False, name=None):
        pieces = [(int(start), int(end))]
        for bs, be in bl_by_chr.get(str(chrom), ()):
            nxt = []
            for s, e in pieces:
                if be <= s or bs >= e:
                    nxt.append((s, e))
                    continue
                if bs > s:
                    nxt.append((s, int(bs)))
                if be < e:
                    nxt.append((int(be), e))
            pieces = nxt
        rows.extend([str(chrom), s, e, *ids] for s, e in pieces if e > s)
    return pd.DataFrame(rows, columns=["#CHR", "START", "END", *ID_COLS])


# configured segmentation, chr-normalized (seg_id = 4th column, else CHR:START-END)
segments = read_BED(segments_bed, extra_columns=("seg_id",))[
    ["#CHR", "START", "END", "seg_id"]
].copy()
assert len(segments) > 0, f"segment_bed, no rows: {segments_bed}"
assert not segments["seg_id"].duplicated().any(), (
    f"segment_bed, duplicate seg_id: {segments_bed}"
)

# chromosome arms; every segment must lie inside exactly one of them
arms = read_BED(region_bed)[["#CHR", "START", "END", "region_id"]]
segments, na_idx = assign_range_to_range(segments, arms, "region_id", rule="contained")
assert len(na_idx) == 0, (
    f"segment_bed, {len(na_idx)} segment(s) not inside one region_bed arm: "
    f"{segments.iloc[na_idx][['#CHR', 'START', 'END']].head(3).to_dict('records')}"
)
logging.info(
    f"{len(segments)} segments over {segments['region_id'].nunique()} arms "
    f"({os.path.basename(segments_bed)})"
)

if blacklist_bed:
    bl = read_BED(blacklist_bed, extra_columns=())[["#CHR", "START", "END"]]
    n_raw = len(segments)
    segments = subtract_blacklist(segments, bl)
    logging.info(
        f"blacklist subtracted: {n_raw} -> {len(segments)} pieces "
        f"({os.path.basename(blacklist_bed)})"
    )

out = segments[["#CHR", "START", "END", *ID_COLS]]
out = out.sort_values(["#CHR", "START"]).reset_index(drop=True)

seg_len = out["END"] - out["START"]
logging.info(
    f"segment bed: {len(out)} rows, {out['seg_id'].nunique()} seg_ids, "
    f"{out['region_id'].nunique()} region_ids; "
    f"length bp: min={int(seg_len.min())}, median={int(seg_len.median())}, "
    f"mean={seg_len.mean():.0f}, max={int(seg_len.max())}"
)

out.to_csv(out_bed, sep="\t", header=False, index=False)
