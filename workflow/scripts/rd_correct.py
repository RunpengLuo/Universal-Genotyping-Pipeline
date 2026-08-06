"""Per-window LOWESS bias correction for bulk samples.

Loads mosdepth fixed-window depth, joins with pre-filtered window BED
(GC/MAP/REPLI/region_id), applies correct_readcount_lowess() per
sample, and saves corrected depth matrix + filtered window dataframe.

The window BED is expected to be pre-filtered by region and blacklist
(produced by build_window_bed.py), with region_id column already present.
"""

import os
import logging

snakemake_handle = snakemake

from utils import (
    set_omp_threads,
    setup_logging,
    add_chr_prefix,
    maybe_path,
    sort_df_chr,
)

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd

from rd_correct_utils import (
    compute_depth_statistics,
    compute_gc_rd_stats,
    correct_readcount_lowess,
    correct_readcount_quadreg,
)
from plot_count_reads import plot_rd_1d_scatter, plot_rd_2d_kde

import matplotlib

matplotlib.use("Agg")
from matplotlib.backends.backend_pdf import PdfPages


# inputs
window_bed_file = snakemake_handle.input["window_bed"]
genome_size = snakemake_handle.input["genome_size"]
region_bed = snakemake_handle.input["region_bed"] or None
blacklist_bed = maybe_path(snakemake_handle.input.get("blacklist_bed", None))

# parameters
assay_type = str(snakemake_handle.params["assay_type"])
sample_id = str(snakemake_handle.params["sample_id"])
dataset_ids = [str(r) for r in snakemake_handle.params["dataset_ids"]]
mosdepth_dir = snakemake_handle.params["mosdepth_dir"]
chroms = snakemake_handle.params["chroms"]
samplesize = int(snakemake_handle.params["samplesize"])
routlier = float(snakemake_handle.params["routlier"])
doutlier = float(snakemake_handle.params["doutlier"])
min_mappability = float(snakemake_handle.params["min_mappability"])
gc_correct = bool(snakemake_handle.params["gc_correct"])
gc_correct_method = str(snakemake_handle.params.get("gc_correct_method", "median"))
rt_correct = bool(snakemake_handle.params["rt_correct"])
qc_dir = snakemake_handle.params["qc_dir"]

# outputs
out_depth_stats = snakemake_handle.output["depth_stats"]
out_dp_corrected = snakemake_handle.output["dp_corrected"]
window_df = snakemake_handle.output["window_df"]

run_id = getattr(snakemake_handle.params, "run_id", "")

nsamples = len(dataset_ids)
sample_ids = [f"{sample_id}_{r}" for r in dataset_ids]
target_chroms = set(chroms)
join_keys = ["#CHR", "START", "END"]

logging.info(f"rd_correct: {nsamples} samples, {len(target_chroms)} chroms")

logging.info("load window BED and mosdepth depth")
win_df = pd.read_table(window_bed_file, sep="\t", dtype={"#CHR": str})
assert "#CHR" in win_df.columns and "GC" in win_df.columns, (
    f"window_bed must have #CHR, START, END, GC columns; got {win_df.columns.tolist()}"
)
win_df["#CHR"] = add_chr_prefix(win_df["#CHR"])

win_df = win_df[win_df["#CHR"].isin(target_chroms)].reset_index(drop=True)

mos_dfs = []
for dataset_id in dataset_ids:
    mos_file = os.path.join(mosdepth_dir, f"{dataset_id}.regions.bed.gz")
    mos_df = pd.read_table(
        mos_file,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "DEPTH"],
        dtype={"#CHR": str},
    )
    mos_df["#CHR"] = add_chr_prefix(mos_df["#CHR"])
    mos_df = mos_df[mos_df["#CHR"].isin(target_chroms)].reset_index(drop=True)
    mos_dfs.append(mos_df)

coords = mos_dfs[0][join_keys].copy()
n_windows = len(coords)
logging.info(f"{n_windows} windows across {len(target_chroms)} chromosomes")

win_df = pd.merge(left=coords, right=win_df, on=join_keys, how="left", sort=False)
_gc_matched = int(win_df["GC"].notna().sum())
logging.info(
    f"GC BED matched {_gc_matched}/{n_windows} ({_gc_matched / max(n_windows, 1) * 100:.1f}%)"
)

dp_raw = np.zeros((n_windows, nsamples), dtype=np.float32)
for i, mos_df in enumerate(mos_dfs):
    dp_raw[:, i] = mos_df["DEPTH"].to_numpy(dtype=np.float32)

# mosdepth emits windows in BAM @SQ order; reorder to genomic order (permute dp_raw too)
win_df["_ord"] = np.arange(len(win_df))
win_df = sort_df_chr(win_df, ch="#CHR", pos="START")
dp_raw = dp_raw[win_df["_ord"].to_numpy()]
win_df = win_df.drop(columns="_ord")

depth_stats = compute_depth_statistics(dp_raw, win_df, sample_ids)
depth_stats.to_csv(out_depth_stats, sep="\t", index=False)
logging.info(f"wrote depth statistics to {out_depth_stats}")
for _, row in depth_stats[depth_stats["#CHR"] == "TOTAL"].iterrows():
    logging.info(
        f"  {row['SAMPLE']}: mean={row['mean_depth']:.2f}, median={row['median_depth']:.2f}"
    )

gc_vals = win_df["GC"].to_numpy()

rd_raw_ylim = max(np.nanquantile(dp_raw, 0.99), 1.0) * 1.1
gc_corr_before, gc_std_before = compute_gc_rd_stats(dp_raw, gc_vals, dataset_ids)

logging.info(f"{n_windows} windows for bias correction")

map_vals = win_df["MAP"].to_numpy() if gc_correct and "MAP" in win_df.columns else None
repli_vals = (
    win_df["REPLI"].to_numpy(dtype=np.float64)
    if rt_correct and "REPLI" in win_df.columns
    else None
)
if repli_vals is not None:
    _n_finite = int(np.isfinite(repli_vals).sum())
    logging.info(
        f"REPLI column: {_n_finite}/{n_windows} ({_n_finite / max(n_windows, 1) * 100:.1f}%) finite"
    )
else:
    logging.info("no REPLI column; skipping replication timing correction")

gc_rmse_list = None
if gc_correct:
    dp_corrected = np.zeros_like(dp_raw, dtype=np.float32)
    gc_rmse_list = []

    if gc_correct_method == "median":
        logging.info("applying correct_readcount_quadreg per sample")
        for i, dataset_id in enumerate(dataset_ids):
            logging.info(f"correcting {dataset_id}")
            dp_corrected[:, i], gc_rmse = correct_readcount_quadreg(
                dp_raw[:, i],
                gc_vals,
                mappability=map_vals,
                repliseq=repli_vals,
                doutlier=doutlier,
                min_mappability=min_mappability,
            )
            gc_rmse_list.append(gc_rmse)
    else:
        logging.info("applying correct_readcount_lowess per sample")
        for i, dataset_id in enumerate(dataset_ids):
            logging.info(f"correcting {dataset_id}")
            dp_corrected[:, i], gc_rmse = correct_readcount_lowess(
                dp_raw[:, i],
                gc_vals,
                mappability=map_vals,
                repliseq=repli_vals,
                samplesize=samplesize,
                routlier=routlier,
                doutlier=doutlier,
                min_mappability=min_mappability,
            )
            gc_rmse_list.append(gc_rmse)
else:
    logging.info("gc_correct=False; skipping bias correction")
    dp_corrected = dp_raw.copy()

rd_ylim = max(np.nanquantile(dp_corrected, 0.99), 1.0) * 1.1

rd_pdf = PdfPages(snakemake_handle.output["qc_pdf"])
plot_rd_1d_scatter(
    win_df,
    dp_raw,
    dp_corrected,
    sample_ids,
    genome_size,
    rd_pdf,
    ylim_before=rd_raw_ylim,
    ylim_after=rd_ylim,
    region_bed=region_bed,
    blacklist_bed=blacklist_bed,
)
plot_rd_2d_kde(
    gc_vals,
    dp_raw,
    dp_corrected,
    sample_ids,
    rd_pdf,
    gc_rmse=gc_rmse_list,
    mappability=map_vals,
    repliseq=repli_vals,
)
rd_pdf.close()

nan_mask = np.isnan(dp_corrected).any(axis=1)
n_nan_rows = int(nan_mask.sum())
n_valid = n_windows - n_nan_rows
logging.info(
    f"NaN row filter: {n_nan_rows}/{n_windows} windows have NaN, "
    f"keeping {n_valid} ({n_valid / max(n_windows, 1) * 100:.1f}%)"
)

if n_nan_rows > 0:
    valid = ~nan_mask
    dp_corrected = dp_corrected[valid]
    win_df = win_df.loc[valid].reset_index(drop=True)

np.savez_compressed(out_dp_corrected, mat=dp_corrected)

out_cols = ["#CHR", "START", "END", "region_id"]
if "seg_id" in win_df.columns:
    out_cols.append("seg_id")
out_cols.append("GC")
if "MAP" in win_df.columns:
    out_cols.append("MAP")
if "REPLI" in win_df.columns:
    out_cols.append("REPLI")
win_df[out_cols].to_csv(
    window_df,
    sep="\t",
    header=True,
    index=False,
    compression="gzip",
)

logging.info("finished rd_correct.")
