"""Build one stream's window BED (bulk), in one pass.

Tiles segment_bed (WGS) or wes_targets (WES), assigns region_id (arm) + seg_id
(breakpoint chunk) by window midpoint, then annotates GC (always), MAP (when a
mappability_bed input is given), and REPLI (when Repli-seq bedGraphs are given).
Output columns: #CHR START END region_id seg_id GC [MAP] [REPLI] -- the stream's
window BED consumed by count_reads. The Repli-seq bigWig fetch + bigWigToBedGraph +
liftOver are Snakemake rules; this script only bins the resulting bedGraphs. Only the
WGS/WES window generators (and the _tile_region they share) are functions.
"""

import logging
import os

snakemake_handle = snakemake  # noqa: F821

t = int(getattr(snakemake_handle, "threads", 1))
for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[var] = str(t)

import numpy as np
import pandas as pd
from pybedtools import BedTool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import get_standard_chroms, read_BED
from plot_utils import _hist_with_stats
from utils import is_canonical_chrom, maybe_path, setup_logging, sort_df_chr

setup_logging(snakemake_handle.log[0])

inp = snakemake_handle.input
p = snakemake_handle.params
out_bed = snakemake_handle.output["window_bed"]
qc_pdf = snakemake_handle.output["qc_pdf"]


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


def generate_wgs_windows(window_size, standard_chroms, region_bed):
    """Tile fixed-size windows within the segment BED (already blacklist-subtracted)."""
    reg_df = pd.read_csv(
        region_bed,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"],
        dtype={"Chromosome": str},
    )
    reg_df = reg_df[reg_df["Chromosome"].isin(standard_chroms)].reset_index(drop=True)
    rows = []
    for _, r in reg_df.iterrows():
        rows.extend(_tile_region(r["Chromosome"], r["Start"], r["End"], window_size))
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


def generate_wes_windows(wes_targets_beds, window_size, standard_chroms):
    """Concatenate the datasets' WES capture targets and adaptively tile them.

    Overlapping targets (within or across the input BEDs) are merged below.
    Blacklist / off-segment windows are dropped by the region_id assignment (their
    midpoint falls outside the segment BED). Genes are not tracked here; they are
    annotated from the GTF (``feature_id``) in combine_counts.
    """
    targets = pd.concat(
        [
            pd.read_csv(f, sep="\t", header=None, comment="#").iloc[:, :3]
            for f in wes_targets_beds
        ],
        ignore_index=True,
    )
    targets.columns = ["Chromosome", "Start", "End"]
    targets = targets[targets["Chromosome"].isin(standard_chroms)].reset_index(
        drop=True
    )
    targets = targets.sort_values(["Chromosome", "Start"]).reset_index(drop=True)

    merged = []
    for chrom, grp in targets.groupby("Chromosome", sort=False):
        starts, ends = grp["Start"].values, grp["End"].values
        cur_start, cur_end = starts[0], ends[0]
        for i in range(1, len(starts)):
            if starts[i] <= cur_end:
                cur_end = max(cur_end, ends[i])
            else:
                merged.append([chrom, cur_start, cur_end])
                cur_start, cur_end = starts[i], ends[i]
        merged.append([chrom, cur_start, cur_end])

    rows = []
    for chrom, start, end in merged:
        rows.extend(_tile_region(chrom, int(start), int(end), window_size))
    return pd.DataFrame(rows, columns=["#CHR", "START", "END"])


standard = get_standard_chroms(p["reference_version"], inp["genome_size"])
region_bed = inp["region_bed"]
genome_size = inp["genome_size"]
logging.info(
    f"build_window_bed: stream={p['mode']}, window_size={p['window_size']}, "
    f"{len(standard)} standard chroms"
)

# tile the segment BED (wgs) or the WES capture targets (wes)
if p["mode"] == "wes":
    windows = generate_wes_windows(
        list(inp["wes_targets_bed"]), int(p["window_size"]), standard
    )
else:
    windows = generate_wgs_windows(int(p["window_size"]), standard, region_bed)
n_tiled = len(windows)
logging.info(f"tiled {n_tiled} windows")

# region_id (arm) + seg_id (chunk) per window by midpoint; drop windows off-segment
regions = read_BED(region_bed)
regions["#CHR"] = regions["#CHR"].astype(str)
regions[["START", "END"]] = regions[["START", "END"]].astype(np.int64)
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

# GC (always): per-window GC fraction via pybedtools nucleotide_content
bt = BedTool.from_dataframe(windows[["#CHR", "START", "END"]])
nuc = bt.nucleotide_content(fi=inp["reference"]).to_dataframe(disable_auto_names=True)
windows["GC"] = nuc["5_pct_gc"].values
logging.info("annotated GC")

# MAP (optional): per-window mean mappability via pybedtools map
mappability_bed = maybe_path(inp["mappability_bed"])
if mappability_bed:
    n_windows = len(windows)
    win_bed = windows[["#CHR", "START", "END"]].copy()
    win_bed["_idx"] = np.arange(n_windows)
    wb = BedTool.from_dataframe(win_bed).sort(g=genome_size)
    map_bt = BedTool(mappability_bed).sort(g=genome_size)
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
        bg = bg[bg["chrom"].isin(standard)].reset_index(drop=True)
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
