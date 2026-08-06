"""Build the bulk window BED, in one pass.

Tiles segment_bed, assigns region_id (arm) + seg_id (breakpoint chunk) by window
midpoint, then annotates GC (always), MAP (when a mappability_bed input is given), and
REPLI (when Repli-seq bedGraphs are given). Output columns: #CHR START END region_id
seg_id GC [MAP] [REPLI] -- the window BED consumed by count_reads (one grid for every
bulk assay: WGS/WGS-lr/WES). The Repli-seq bigWig fetch + bigWigToBedGraph + liftOver
are Snakemake rules; this script only bins the resulting bedGraphs.
"""

import logging
import os

snakemake_handle = snakemake  # noqa: F821

from utils import (
    set_omp_threads,
    setup_logging,
    is_canonical_chrom,
    match_chr_style,
    maybe_path,
    sort_df_chr,
    strip_chr_prefix,
)

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from pybedtools import BedTool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import read_BED
from plot_utils import _hist_with_stats


inp = snakemake_handle.input
p = snakemake_handle.params
out_bed = snakemake_handle.output["window_bed"]
qc_pdf = snakemake_handle.output["qc_pdf"]


def _to_genome_style(feature):
    """Rename one pybedtools interval to the genome's naming."""
    feature.chrom = match_chr_style(feature.chrom, input_nochr)
    return feature


def generate_wgs_windows(window_size, chroms, regions):
    """Tile fixed-size windows within the segment BED (already blacklist-subtracted)."""

    def _tile_region(chrom, start, end, window_size):
        """Tile a region into fixed-size windows, merging an undersized last bin."""
        rows = []
        pos = start
        while pos < end:
            w_end = min(pos + window_size, end)
            rows.append([chrom, pos, w_end])
            pos = w_end
        if len(rows) > 1 and (rows[-1][2] - rows[-1][1]) < window_size // 2:
            rows[-2][2] = rows[-1][2]
            rows.pop()
        return rows

    reg_df = regions[regions["#CHR"].isin(chroms)].reset_index(drop=True)
    rows = []
    for _, r in reg_df.iterrows():
        rows.extend(_tile_region(r["#CHR"], r["START"], r["END"], window_size))
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


region_bed = inp["region_bed"]
genome_size = inp["genome_size"]
chroms = list(p["chroms"])
input_nochr = p["input_nochr"]
logging.info(f"build_window_bed: window_size={p['window_size']}, {len(chroms)} chroms")

regions = read_BED(region_bed)
regions["#CHR"] = regions["#CHR"].astype(str)
regions[["START", "END"]] = regions[["START", "END"]].astype(np.int64)

# tile the segment BED (one grid for every bulk assay: WGS/WGS-lr/WES)
windows = generate_wgs_windows(int(p["window_size"]), chroms, regions)
n_tiled = len(windows)
logging.info(f"tiled {n_tiled} windows")

# region_id (arm) + seg_id (chunk) per window by midpoint; drop windows off-segment
mids = (windows["START"] + windows["END"]) // 2
region_out = pd.Series(pd.NA, index=windows.index, dtype="object")
seg_out = pd.Series(pd.NA, index=windows.index, dtype="object")
for chrom in regions["#CHR"].unique():
    win_mask = windows["#CHR"] == chrom
    if not win_mask.any():
        continue
    reg_ch = regions.loc[regions["#CHR"] == chrom].sort_values("START")
    starts = reg_ch["START"].to_numpy()
    ends = reg_ch["END"].to_numpy()
    ids = reg_ch["region_id"].to_numpy()
    segs = reg_ch["seg_id"].to_numpy()
    positions = mids[win_mask].to_numpy()
    idx = np.searchsorted(starts, positions, side="right") - 1
    valid = (idx >= 0) & (positions < ends[idx.clip(min=0)])
    win_indices = windows.index[win_mask]
    region_out.loc[win_indices[valid]] = ids[idx[valid]]
    seg_out.loc[win_indices[valid]] = segs[idx[valid]]
windows["region_id"] = region_out
windows["seg_id"] = seg_out
windows = windows[windows["region_id"].notna()].reset_index(drop=True)
logging.info(
    f"region_id/seg_id: kept {len(windows)}/{n_tiled} windows (dropped off-segment)"
)

# drop non-canonical contigs, then sort genomically
n_pre = len(windows)
windows = windows[windows["#CHR"].map(is_canonical_chrom)].copy()
windows = sort_df_chr(windows, ch="#CHR", pos="START")
logging.info(
    f"sorted {len(windows)} windows on {windows['#CHR'].nunique()} chroms "
    f"(dropped {n_pre - len(windows)} non-canonical)"
)

# bedtools resolves contigs against the reference FASTA and genome_size, so the
# intervals it receives carry the genome's naming; results map back by row order/_idx
bed_windows = windows[["#CHR", "START", "END"]].copy()
if input_nochr:
    bed_windows["#CHR"] = bed_windows["#CHR"].map(strip_chr_prefix)

# GC (always): per-window GC fraction via pybedtools nucleotide_content
bt = BedTool.from_dataframe(bed_windows)
nuc = bt.nucleotide_content(fi=inp["reference"]).to_dataframe(disable_auto_names=True)
windows["GC"] = nuc["5_pct_gc"].values
logging.info("annotated GC")

# MAP (optional): per-window mean mappability via pybedtools map
mappability_bed = maybe_path(inp["mappability_bed"])
if mappability_bed:
    n_windows = len(windows)
    win_bed = bed_windows.copy()
    win_bed["_idx"] = np.arange(n_windows)
    wb = BedTool.from_dataframe(win_bed).sort(g=genome_size)
    # streams; bedtools resolves both operands against genome_size
    map_bt = (
        BedTool(mappability_bed).each(_to_genome_style).saveas().sort(g=genome_size)
    )
    map_cov = pd.read_csv(
        wb.map(b=map_bt, c=4, o="mean", g=genome_size).fn,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "_idx", "MAP"],
        dtype={"#CHR": str, "START": int, "END": int, "_idx": int, "MAP": str},
    )
    assert len(map_cov) == n_windows, (
        f"bedtools map returned {len(map_cov)} rows, expected {n_windows}"
    )
    map_cov["MAP"] = (
        pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    )
    windows["MAP"] = map_cov.sort_values("_idx")["MAP"].values
    logging.info(f"annotated MAP from {os.path.basename(mappability_bed)}")
else:
    logging.info("no mappability_bed; MAP skipped")

# REPLI (optional): bin each Repli-seq bedGraph by midpoint, average across tracks
bedgraphs = list(inp["bedgraphs"])
if bedgraphs:
    n_win = len(windows)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int32)
    for lf in bedgraphs:
        bg = pd.read_csv(
            lf,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "signal"],
            dtype={
                "chrom": str,
                "start": np.int64,
                "end": np.int64,
                "signal": np.float64,
            },
        )
        bg = bg[bg["chrom"].isin(chroms)].reset_index(drop=True)
        for chrom, grp in bg.groupby("chrom", sort=False):
            win_idx = np.where((windows["#CHR"] == chrom).to_numpy())[0]
            if len(win_idx) == 0:
                continue
            win_starts = windows["START"].to_numpy()[win_idx]
            order = np.argsort(win_starts)
            win_starts, win_idx = win_starts[order], win_idx[order]
            bg_mids = ((grp["start"] + grp["end"]) // 2).to_numpy()
            bin_idx = np.searchsorted(win_starts, bg_mids, side="right") - 1
            valid = (bin_idx >= 0) & (bin_idx < len(win_idx))
            global_idx = win_idx[bin_idx[valid]]
            signal_sum[global_idx] += grp["signal"].to_numpy()[valid]
            signal_count[global_idx] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        windows["REPLI"] = np.round(
            np.where(signal_count > 0, signal_sum / signal_count, np.nan), 6
        )
    logging.info(f"annotated REPLI from {len(bedgraphs)} bedGraph track(s)")
else:
    logging.info("no Repli-seq bedGraphs; REPLI skipped")

windows.to_csv(out_bed, sep="\t", index=False, compression="gzip")
logging.info(
    f"wrote {len(windows)} windows [{', '.join(windows.columns)}] -> {out_bed}"
)

# QC: segment-length (segment.bed) and window-length distributions
seg_len_kbp = (regions["END"] - regions["START"]).to_numpy() / 1000.0
win_len = (windows["END"] - windows["START"]).to_numpy()
fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 4))
_hist_with_stats(
    ax0, seg_len_kbp, "segment length (kbp)", header="segment.bed", ylabel="# segments"
)
_hist_with_stats(
    ax1,
    win_len,
    "window length (bp)",
    header=f"{p['mode']} windows",
    ylabel="# windows",
)
fig.tight_layout()
with PdfPages(qc_pdf) as pdf:
    pdf.savefig(fig)
plt.close(fig)
logging.info(f"wrote QC length histograms -> {qc_pdf}")
