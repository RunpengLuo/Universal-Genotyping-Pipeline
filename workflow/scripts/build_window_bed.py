"""Build the window BED, in one pass.

Tiles segment_bed per segment row (so no window spans two segments), assigns region_id
(arm) + seg_id (segment) by window midpoint, then annotates each bias-correction
covariate whose input is non-empty: GC (reference FASTA), MAP (mappability_bed), REPLI
(Repli-seq bedGraphs). Output columns: #CHR START END region_id seg_id [GC] [MAP]
[REPLI] -- one bin set for every assay of the run. Bulk feeds all three covariates and
counts the windows with mosdepth -> rd_correct; single-cell feeds none of them and uses
the intervals as its fixed-bin skeleton. The Repli-seq bigWig fetch + bigWigToBedGraph +
liftOver are Snakemake rules; this script only bins the resulting bedGraphs.
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
from range_utils import assign_range_to_range


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

# tile the segment BED (one bin set for every bulk assay: WGS/WGS-lr/WES)
windows = generate_wgs_windows(int(p["window_size"]), chroms, regions)
n_tiled = len(windows)
logging.info(f"tiled {n_tiled} windows")

# region_id (arm) + seg_id (segment) per window by midpoint; drop windows off-segment
windows, _ = assign_range_to_range(windows, regions, "region_id", rule="midpoint")
windows, _ = assign_range_to_range(
    windows, regions, "seg_id", rule="midpoint", dropna=True
)
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
# intervals it receives carry the genome's naming; results map back by row order/_df_idx
bed_windows = windows[["#CHR", "START", "END"]].copy()
if input_nochr:
    bed_windows["#CHR"] = bed_windows["#CHR"].map(strip_chr_prefix)

# GC (optional): per-window GC fraction via pybedtools nucleotide_content
reference = maybe_path(inp["reference"])
if reference:
    bt = BedTool.from_dataframe(bed_windows)
    nuc = bt.nucleotide_content(fi=reference).to_dataframe(disable_auto_names=True)
    windows["GC"] = nuc["5_pct_gc"].values
    logging.info("annotated GC")
else:
    logging.info("no reference; GC skipped")

# MAP (optional): per-window mean mappability via pybedtools map
mappability_bed = maybe_path(inp["mappability_bed"])
if mappability_bed:
    n_windows = len(windows)
    win_bed = bed_windows.copy()
    win_bed["_df_idx"] = np.arange(n_windows)
    wb = BedTool.from_dataframe(win_bed).sort(g=genome_size)
    # streams; bedtools resolves both operands against genome_size
    map_bt = (
        BedTool(mappability_bed).each(_to_genome_style).saveas().sort(g=genome_size)
    )
    map_cov = pd.read_csv(
        wb.map(b=map_bt, c=4, o="mean", g=genome_size).fn,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "_df_idx", "MAP"],
        dtype={"#CHR": str, "START": int, "END": int, "_df_idx": int, "MAP": str},
    )
    assert len(map_cov) == n_windows, (
        f"bedtools map, got {len(map_cov)} rows, expected {n_windows}"
    )
    map_cov["MAP"] = (
        pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    )
    windows["MAP"] = map_cov.sort_values("_df_idx")["MAP"].values
    logging.info(f"annotated MAP from {os.path.basename(mappability_bed)}")
else:
    logging.info("no mappability_bed; MAP skipped")

# REPLI (optional): bin each Repli-seq bedGraph by midpoint, average across tracks
bedgraphs = list(inp["bedgraphs"])
if bedgraphs:
    n_win = len(windows)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int64)
    win_ranges = windows[["#CHR", "START", "END"]].assign(_win=np.arange(n_win))
    for lf in bedgraphs:
        bg = pd.read_csv(
            lf,
            sep="\t",
            header=None,
            names=["#CHR", "START", "END", "signal"],
            dtype={
                "#CHR": str,
                "START": np.int64,
                "END": np.int64,
                "signal": np.float64,
            },
        )
        bg = bg[bg["#CHR"].isin(chroms)].reset_index(drop=True)
        bg, _ = assign_range_to_range(
            bg, win_ranges, "_win", rule="midpoint", dropna=True
        )
        # bincount, not `arr[idx] +=`: several bedGraph rows can share one window
        hit = bg["_win"].to_numpy()
        signal_sum += np.bincount(hit, weights=bg["signal"].to_numpy(), minlength=n_win)
        signal_count += np.bincount(hit, minlength=n_win)
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
    header=f"{p['window_size']} bp windows",
    ylabel="# windows",
)
fig.tight_layout()
with PdfPages(qc_pdf) as pdf:
    pdf.savefig(fig)
plt.close(fig)
logging.info(f"wrote QC length histograms -> {qc_pdf}")
