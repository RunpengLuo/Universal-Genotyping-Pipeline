"""SNP-informed adaptive binning + window→bin depth aggregation + RDR.

Combines:
1. Adaptive binning from SNP data (switch probs, MSR/MSPB)
2. Aggregation of corrected window depth into adaptive bins
3. RDR computation (tumor/normal ratio or median-centering)
"""

import os
import shutil
import logging

t = int(getattr(snakemake, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import setup_logging, maybe_path, stamp_path
from aggregation_utils import (
    adaptive_binning_windows,
    assign_pos_to_range,
    detect_phase_flips,
    matrix_segmentation,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_allele_freqs, plot_rdr_baf
from switchprobs import (
    interp_cM_blocks,
    estimate_switchprobs_cM,
    estimate_switchprobs_PS,
)

setup_logging(snakemake.log[0])

snp_info = snakemake.input["snp_info"]
tot_mtx_snp = snakemake.input["tot_mtx_snp"]
a_mtx_snp = snakemake.input["a_mtx_snp"]
b_mtx_snp = snakemake.input["b_mtx_snp"]
dp_corrected = snakemake.input["dp_corrected"]
window_df_file = snakemake.input["window_df"]

gmap_file = maybe_path(snakemake.input["gmap_file"])
region_bed = snakemake.input["region_bed"]
blacklist_bed = maybe_path(snakemake.input.get("blacklist_bed", None))
genome_size = snakemake.input["genome_size"]
gtf_file = maybe_path(snakemake.input["gtf_file"])

qc_dir = snakemake.params["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(snakemake.params, "run_id", "")

sample_df = pd.read_table(snakemake.input["sample_file"])
sample_name = sample_df["SAMPLE"].iloc[0]
rep_ids = sample_df["REP_ID"].tolist()
sample_types = sample_df["sample_type"].tolist()
assay_type = snakemake.params["assay_type"]
chromosomes = snakemake.params["chromosomes"]

has_normal = "normal" in sample_types
tumor_sidx = 1 if has_normal else 0
nsamples = len(sample_df)

logging.info(
    f"combine_counts: sample={sample_name}, assay={assay_type}, "
    f"nsamples={nsamples}, has_normal={has_normal}"
)

snps = pd.read_table(snp_info, sep="\t")
tot_mtx = np.load(tot_mtx_snp)["mat"].astype(np.int32)
a_mtx = np.load(a_mtx_snp)["mat"].astype(np.int32)
b_mtx = np.load(b_mtx_snp)["mat"].astype(np.int32)
nsnps = tot_mtx.shape[0]
logging.info(f"loaded {nsnps} SNPs, {nsamples} samples")

dp_mtx = np.load(dp_corrected)["mat"]
window_df = pd.read_table(window_df_file, sep="\t")
n_window_df = len(window_df)
logging.info(f"loaded {n_window_df} corrected window_df")

assert "region_id" in snps.columns, "invalid SNP file"
grp_cols = ["region_id"]
if "PS" not in snps.columns:
    logging.info("PS not in SNP columns, setting PS=1 for all SNPs")
    snps["PS"] = 1
else:
    logging.info("PS is provided")
num_phaseset = snps["PS"].nunique()
logging.info(f"#phaseset={num_phaseset}")
grp_cols.append("PS")

phase_flip_test = bool(snakemake.params["phase_flip_test"])
if phase_flip_test:
    snps["phase_group"] = detect_phase_flips(
        snps,
        a_mtx,
        b_mtx,
        grp_cols=grp_cols,
        tumor_sidx=tumor_sidx,
        epsilon=float(snakemake.params["phase_flip_epsilon"]),
        alpha=float(snakemake.params["phase_flip_alpha"]),
    )
    grp_cols.append("phase_group")

window_df["win_idx"] = np.arange(len(window_df))
# Assign PS (and phase_group if active) to window_df from SNPs via majority vote.
# Only the columns needed for window annotation are extracted to avoid copying
# the full SNP DataFrame.
snp_window_cols = ["#CHR", "POS0", "PS"]
if phase_flip_test:
    snp_window_cols.append("phase_group")
_snps_tmp = snps[snp_window_cols].copy()
_snps_tmp = assign_pos_to_range(_snps_tmp, window_df, ref_id="win_idx", pos_col="POS0")
_snps_tmp = _snps_tmp.dropna(subset=["win_idx"])
_snps_tmp["win_idx"] = _snps_tmp["win_idx"].astype(np.int64)
snps_per_win = _snps_tmp.groupby("win_idx").size()
n_win_with_snps = len(snps_per_win)
n_win_total = len(window_df)
logging.info(
    f"SNPs per window: {n_win_with_snps}/{n_win_total} windows have SNPs, "
    f"mean={snps_per_win.mean():.1f}, median={snps_per_win.median():.1f}, "
    f"min={snps_per_win.min()}, max={snps_per_win.max()}"
)
win_ps = _snps_tmp.groupby("win_idx")["PS"].agg(lambda x: x.mode().iloc[0])
window_df["PS"] = window_df["win_idx"].map(win_ps)
# Windows with no SNPs inherit PS from previous window (forward-fill)
if window_df["PS"].isna().any():
    window_df["PS"] = window_df["PS"].ffill()

if phase_flip_test:
    win_pg = _snps_tmp.groupby("win_idx")["phase_group"].agg(lambda x: x.mode().iloc[0])
    window_df["phase_group"] = window_df["win_idx"].map(win_pg)
    if window_df["phase_group"].isna().any():
        window_df["phase_group"] = window_df["phase_group"].ffill()

max_blocksize = int(snakemake.params["max_blocksize"])
bbs, snps = adaptive_binning_windows(
    window_df,
    snps,
    tot_mtx,
    int(snakemake.params["min_snp_reads"]),
    int(snakemake.params["min_snp_per_block"]),
    grp_cols=grp_cols,
    tumor_sidx=tumor_sidx,
    max_blocksize=max_blocksize,
)
num_bbs = len(bbs)
bb_ids = snps["bb_id"].to_numpy()

snp_orig_idx = snps["_orig_idx"].to_numpy()
tot_mtx = tot_mtx[snp_orig_idx]
a_mtx = a_mtx[snp_orig_idx]
b_mtx = b_mtx[snp_orig_idx]

a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, num_bbs)
b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, num_bbs)
tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, num_bbs)

pdf_path = stamp_path(os.path.join(qc_dir, "combine_counts.pdf"), run_id)

baf_mtx_bb = np.divide(
    b_mtx_bb,
    tot_mtx_bb,
    where=tot_mtx_bb > 0,
    out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
)

logging.info("aggregating corrected window depth into adaptive bins")

win_bin_ids = window_df["bin_id"].to_numpy()
win_lengths = (window_df["END"] - window_df["START"]).to_numpy(dtype=np.float64)
total_len_per_bin = np.bincount(win_bin_ids, weights=win_lengths, minlength=num_bbs)

bb_dp = np.full((num_bbs, nsamples), np.nan, dtype=np.float32)
for s in range(nsamples):
    weighted_sums = np.bincount(
        win_bin_ids, weights=dp_mtx[:, s] * win_lengths, minlength=num_bbs
    )
    with np.errstate(invalid="ignore"):
        bb_dp[:, s] = weighted_sums / total_len_per_bin

median_normalization = bool(getattr(snakemake.params, "median_normalization", False))
use_normal = has_normal and not median_normalization
logging.info(f"compute bb RDR, use_normal={use_normal}")

if use_normal:
    window_sizes = (window_df["END"] - window_df["START"]).to_numpy(dtype=np.float64)
    total_bases = np.nansum(dp_mtx * window_sizes[:, None], axis=0)
    library_correction = total_bases[0] / total_bases[1:]
    logging.info(f"library normalization factor: {library_correction}")

    normal_bb_dp = bb_dp[:, 0]
    with np.errstate(invalid="ignore", divide="ignore"):
        bb_rdr = bb_dp[:, tumor_sidx:] / normal_bb_dp[:, None]
        bb_rdr *= library_correction[None, :]
else:
    rdr_dp = bb_dp[:, tumor_sidx:]
    rdr_reps = rep_ids[tumor_sidx:]
    n_rdr = rdr_dp.shape[1]
    bb_rdr = np.full((num_bbs, n_rdr), np.nan, dtype=np.float32)
    for i in range(n_rdr):
        col = rdr_dp[:, i]
        valid_i = np.isfinite(col) & (col > 0)
        if valid_i.any():
            med = np.median(col[valid_i])
            logging.info(f"  bb median-centering {rdr_reps[i]}: median={med:.4f}")
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[valid_i, i] = col[valid_i] / med

rdr_outlier_quantile = float(snakemake.params["rdr_outlier_quantile"])
if rdr_outlier_quantile > 0:
    rdr_upper = np.nanquantile(bb_rdr, 1 - rdr_outlier_quantile)
    n_outlier = int(np.nansum(bb_rdr > rdr_upper))
    logging.info(
        f"RDR outlier filter: quantile={rdr_outlier_quantile}, "
        f"threshold={rdr_upper:.4f}, {n_outlier} entries set to NaN"
    )
    bb_rdr[bb_rdr > rdr_upper] = np.nan

n_nan_bb = int(np.isnan(bb_rdr[:, 0]).sum())
_n_filled = num_bbs - n_nan_bb
logging.info(
    f"bb RDR: {_n_filled}/{num_bbs} ({_n_filled / max(num_bbs, 1) * 100:.1f}%) filled, "
    f"{n_nan_bb} NaN"
)

tumor_rep_ids = rep_ids[tumor_sidx:]
rdr_ylim = (np.round(np.nanquantile(bb_rdr, 0.99)).astype(int) + 1) * 1.1

with PdfPages(pdf_path) as pdf:
    plot_allele_freqs(
        bbs,
        rep_ids,
        tot_mtx_bb,
        b_mtx_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=False,
        allele="B",
        unit="bb",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        pdf=pdf,
    )
    plot_rdr_baf(
        bbs,
        bb_rdr,
        baf_mtx_bb[:, tumor_sidx:],
        list(tumor_rep_ids),
        genome_size,
        pdf_path,
        unit="bb",
        rdr_ylim=rdr_ylim,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        pdf=pdf,
    )
logging.info(f"saved QC PDF to {pdf_path}")

nan_mask = (
    np.isnan(baf_mtx_bb).any(axis=1)
    | np.isnan(bb_dp).any(axis=1)
    | np.isnan(bb_rdr).any(axis=1)
)
n_nan_rows = int(nan_mask.sum())
n_valid = num_bbs - n_nan_rows
logging.info(
    f"NaN row filter: {n_nan_rows}/{num_bbs} bins have NaN, "
    f"keeping {n_valid} ({n_valid / max(num_bbs, 1) * 100:.1f}%)"
)

if n_nan_rows > 0:
    valid = ~nan_mask
    bbs = bbs.loc[valid].reset_index(drop=True)
    tot_mtx_bb = tot_mtx_bb[valid]
    a_mtx_bb = a_mtx_bb[valid]
    b_mtx_bb = b_mtx_bb[valid]
    baf_mtx_bb = baf_mtx_bb[valid]
    bb_dp = bb_dp[valid]
    bb_rdr = bb_rdr[valid]

kept_bb_ids = np.where(~nan_mask)[0] if n_nan_rows > 0 else np.arange(num_bbs)
old_to_new = {old: new for new, old in enumerate(kept_bb_ids)}
snps_valid = snps[snps["bb_id"].isin(old_to_new)].copy()
snps_valid["bb_id"] = snps_valid["bb_id"].map(old_to_new)
bbs["bb_id"] = np.arange(len(bbs))

logging.info("estimate bin-level switchprobs")
if gmap_file is not None:
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms = interp_cM_blocks(bbs, snps_valid, genetic_map, block_id_col="bb_id")
    bbs["switchprobs"] = estimate_switchprobs_cM(
        dist_cms,
        nu=float(snakemake.params["nu"]),
        min_switchprob=float(snakemake.params["min_switchprob"]),
    )
else:
    switchprob_ps = float(snakemake.params["switchprob_ps"])
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bbs[["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]].to_csv(
    snakemake.output["bb_file"], sep="\t", header=True, index=False
)
np.savez_compressed(snakemake.output["tot_mtx_bb"], mat=tot_mtx_bb)
np.savez_compressed(snakemake.output["a_mtx_bb"], mat=a_mtx_bb)
np.savez_compressed(snakemake.output["b_mtx_bb"], mat=b_mtx_bb)
np.savez_compressed(snakemake.output["baf_mtx_bb"], mat=baf_mtx_bb)
np.savez_compressed(snakemake.output["dp_mtx_bb"], mat=bb_dp)
np.savez_compressed(snakemake.output["rdr_mtx_bb"], mat=bb_rdr)
shutil.copy2(snakemake.input["sample_file"], snakemake.output["sample_file"])
logging.info("finished combine_counts.")
