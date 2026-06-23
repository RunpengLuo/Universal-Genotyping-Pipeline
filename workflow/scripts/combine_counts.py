"""SNP-informed adaptive binning across all bulk assays + depth aggregation + RDR.

All bulk assays present in the run (e.g. bulkWGS + bulkWGS-lr) are segmented on ONE
shared bin grid: ``adaptive_segmentation`` requires ``min_snp_reads`` in every tumor
sample, so stacking all assays' tumor columns yields bins that jointly satisfy every
bulk sample. Allele counts are aggregated per bin across all samples; read depth and
RDR are computed per assay, normalizing each assay's tumors by that assay's own normal.

The allele matrices come as one joint set from phase_and_concat_bulk (read directly, no
union); depth/window inputs stay per-assay (index-aligned to ``params.bulk_assays``).
Outputs go under ``bb_dir/{stream}/`` (stream = bulkWGS or bulkWES); matrix columns are
the bulk samples.
"""

import os
import logging

snakemake_handle = snakemake

t = int(getattr(snakemake_handle, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from utils import setup_logging, maybe_path, sort_df_chr
from aggregation_utils import (
    adaptive_segmentation,
    assign_pos_to_range,
    count_split_genes,
    detect_phase_flips,
    gene_block_labels,
    matrix_segmentation,
    merge_feature_ids,
)
from combine_counts_utils import (
    aggregate_window_depth_to_bins,
    build_assay_blocks,
    compute_bb_rdr,
    load_bulk_snp_matrices,
    setup_phaseset_groups,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_allele_freqs, plot_rdr_baf, plot_segmentation_qc
from switchprobs import (
    interp_cM_blocks,
    estimate_switchprobs_cM,
    estimate_switchprobs_PS,
)

log_file = snakemake_handle.log[0]
setup_logging(log_file)

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
dp_corrected_files = list(snakemake_handle.input["dp_corrected"])
window_df_files = list(snakemake_handle.input["window_df"])
sample_file = snakemake_handle.input["sample_file"]
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = maybe_path(snakemake_handle.input.get("blacklist_bed", None))
genome_size = snakemake_handle.input["genome_size"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
run_id = snakemake_handle.params["run_id"]
bulk_assays = list(snakemake_handle.params["bulk_assays"])
median_normalization = bool(snakemake_handle.params["median_normalization"])
phase_flip_test = bool(snakemake_handle.params["phase_flip_test"])
phase_flip_epsilon = float(snakemake_handle.params["phase_flip_epsilon"])
phase_flip_alpha = float(snakemake_handle.params["phase_flip_alpha"])
gene_aware_binning_param = bool(snakemake_handle.params["gene_aware_binning"])
max_blocksize = int(snakemake_handle.params["max_blocksize"])
min_snp_reads = int(snakemake_handle.params["min_snp_reads"])
min_snp_per_block = int(snakemake_handle.params["min_snp_per_block"])
rdr_outlier_quantile = float(snakemake_handle.params["rdr_outlier_quantile"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])

# outputs
out_bb_file = snakemake_handle.output["bb_file"]
out_tot_mtx_bb = snakemake_handle.output["tot_mtx_bb"]
out_a_mtx_bb = snakemake_handle.output["a_mtx_bb"]
out_b_mtx_bb = snakemake_handle.output["b_mtx_bb"]
out_dp_mtx_bb = snakemake_handle.output["dp_mtx_bb"]
out_rdr_mtx_bb = snakemake_handle.output["rdr_mtx_bb"]
out_sample_file = snakemake_handle.output["sample_file"]

sample_df = pd.read_table(sample_file)
snps, tot_mtx, a_mtx, b_mtx = load_bulk_snp_matrices(
    snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp
)
dp_corrected_list = [np.load(f)["mat"] for f in dp_corrected_files]
window_df_list = [pd.read_table(f, sep="\t") for f in window_df_files]
n_snps = len(snps)

sample_name = sample_df["SAMPLE_NAME"].iloc[0]
total_samples = len(sample_df)
logging.info(
    f"combine_counts: sample={sample_name}, bulk_assays={bulk_assays}; "
    f"{n_snps} SNPs x {total_samples} samples"
)

col_assay = sample_df["assay_type"].tolist()
col_repid = sample_df["REP_ID"].tolist()
assay_blocks, tumor_cols_all = build_assay_blocks(sample_df, bulk_assays)
logging.info(f"{total_samples} bulk samples, {len(tumor_cols_all)} tumor columns")

has_feature = "feature_id" in snps.columns
grp_cols = setup_phaseset_groups(snps)

if phase_flip_test:
    snps["phase_group"] = detect_phase_flips(
        snps,
        a_mtx[:, tumor_cols_all],
        b_mtx[:, tumor_cols_all],
        grp_cols=grp_cols,
        tumor_sidx=0,
        epsilon=phase_flip_epsilon,
        alpha=phase_flip_alpha,
    )
    grp_cols.append("phase_group")

window_df = (
    pd.concat([w[["#CHR", "START", "END", "region_id"]] for w in window_df_list], ignore_index=True)
    .drop_duplicates(["#CHR", "START", "END"])
)
window_df = sort_df_chr(window_df, ch="#CHR", pos="START").reset_index(drop=True)

gene_aware_binning = gene_aware_binning_param and has_feature
window_df["win_idx"] = np.arange(len(window_df))
snp_window_cols = ["#CHR", "POS0", "PS"]
if phase_flip_test:
    snp_window_cols.append("phase_group")
if gene_aware_binning:
    snp_window_cols.append("feature_id")
_snps_tmp = snps[snp_window_cols].copy()
_snps_tmp = assign_pos_to_range(_snps_tmp, window_df, ref_id="win_idx", pos_col="POS0")
_snps_tmp = _snps_tmp.dropna(subset=["win_idx"])
_snps_tmp["win_idx"] = _snps_tmp["win_idx"].astype(np.int64)
snps_per_win = _snps_tmp.groupby("win_idx").size()
logging.info(
    f"SNPs per window: {len(snps_per_win)}/{len(window_df)} windows have SNPs, "
    f"mean={snps_per_win.mean():.1f}, median={snps_per_win.median():.1f}"
)
win_ps = _snps_tmp.groupby("win_idx")["PS"].agg(lambda x: x.mode().iloc[0])
window_df["PS"] = window_df["win_idx"].map(win_ps)
if window_df["PS"].isna().any():
    window_df["PS"] = window_df["PS"].ffill()

if phase_flip_test:
    win_pg = _snps_tmp.groupby("win_idx")["phase_group"].agg(lambda x: x.mode().iloc[0])
    window_df["phase_group"] = window_df["win_idx"].map(win_pg)
    if window_df["phase_group"].isna().any():
        window_df["phase_group"] = window_df["phase_group"].ffill()

if gene_aware_binning:
    # glue each gene's window span into one block so a bin never splits a gene;
    # explode the ;-joined multi-gene feature_id so each gene gets its own span
    genic = _snps_tmp[
        _snps_tmp["feature_id"].notna() & (_snps_tmp["feature_id"] != "intergenic")
    ].copy()
    genic["feature_id"] = genic["feature_id"].str.split(";")
    genic = genic.explode("feature_id")
    genic = genic[genic["feature_id"] != "intergenic"]
    rng = genic.groupby("feature_id")["win_idx"].agg(["min", "max"])
    window_df["gene_block"] = gene_block_labels(
        len(window_df), zip(rng["min"].to_numpy(), rng["max"].to_numpy())
    )
    logging.info(
        f"gene-aware binning: {len(rng)} genes over {len(window_df)} windows -> "
        f"{window_df['gene_block'].nunique()} gene/intergenic blocks (bins never split a gene)"
    )

tot_tumor = np.ascontiguousarray(tot_mtx[:, tumor_cols_all], dtype=np.float64)
bbs, snps = adaptive_segmentation(
    window_df,
    snps,
    tot_tumor,
    min_snp_reads,
    min_snp_per_block,
    grp_cols=grp_cols,
    tumor_sidx=0,
    max_blocksize=max_blocksize,
    gene_aware=gene_aware_binning,
)
num_bbs = len(bbs)

count_split_genes(snps, grp_cols, gene_aware_binning)

bb_ids = snps["bb_id"].to_numpy()

snp_orig_idx = snps["_orig_idx"].to_numpy()
tot_mtx = tot_mtx[snp_orig_idx]
a_mtx = a_mtx[snp_orig_idx]
b_mtx = b_mtx[snp_orig_idx]

a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, num_bbs)
b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, num_bbs)
tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, num_bbs)

baf_mtx_bb = np.divide(
    b_mtx_bb,
    tot_mtx_bb,
    where=tot_mtx_bb > 0,
    out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
)

logging.info("aggregating corrected window depth into adaptive bins (per assay)")
bb_dp, bb_bases = aggregate_window_depth_to_bins(
    assay_blocks, window_df, window_df_list, dp_corrected_list, num_bbs, total_samples
)

logging.info(f"compute bb RDR, median_normalization={median_normalization}")
bb_rdr = compute_bb_rdr(
    assay_blocks,
    window_df_list,
    dp_corrected_list,
    bb_dp,
    tumor_cols_all,
    median_normalization,
    rdr_outlier_quantile,
    col_repid,
)

sample_labels = [f"{col_assay[i]}:{col_repid[i]}" for i in range(total_samples)]
tumor_labels = [sample_labels[c] for c in tumor_cols_all]
rdr_ylim = (np.round(np.nanquantile(bb_rdr, 0.99)).astype(int) + 1) * 1.1

if gene_aware_binning and "feature_id" in snps.columns:
    _genic = snps[
        snps["feature_id"].notna() & (snps["feature_id"] != "intergenic")
    ][["bb_id", "feature_id"]].copy()
    _genic["feature_id"] = _genic["feature_id"].str.split(";")
    _genic = _genic.explode("feature_id")
    _genic = _genic[_genic["feature_id"] != "intergenic"]
    bb_gene_count = (
        _genic.groupby("bb_id")["feature_id"].nunique()
        .reindex(range(num_bbs)).fillna(0).to_numpy()
    )
else:
    bb_gene_count = None

pdf_path = snakemake_handle.output["qc_pdf"]
with PdfPages(pdf_path) as pdf:
    plot_segmentation_qc(
        bbs,
        sample_df,
        bb_bases,
        b_mtx_bb,
        tot_mtx_bb,
        pdf=pdf,
        gene_count=bb_gene_count,
    )
    plot_allele_freqs(
        bbs,
        sample_labels,
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
        baf_mtx_bb[:, tumor_cols_all],
        tumor_labels,
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
        nu=nu,
        min_switchprob=min_switchprob,
    )
else:
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bb_cols = ["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]
if "feature_id" in snps_valid.columns:
    bbs["feature_id"] = (
        bbs["bb_id"]
        .map(snps_valid.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )
    bb_cols.append("feature_id")
bbs[bb_cols].to_csv(out_bb_file, sep="\t", header=True, index=False)
np.savez_compressed(out_tot_mtx_bb, mat=tot_mtx_bb)
np.savez_compressed(out_a_mtx_bb, mat=a_mtx_bb)
np.savez_compressed(out_b_mtx_bb, mat=b_mtx_bb)
np.savez_compressed(out_dp_mtx_bb, mat=bb_dp)
np.savez_compressed(out_rdr_mtx_bb, mat=bb_rdr)

sample_df.to_csv(out_sample_file, sep="\t", index=False)
logging.info("finished combine_counts.")
