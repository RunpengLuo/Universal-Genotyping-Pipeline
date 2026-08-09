"""SNP-informed adaptive binning across all bulk assays + depth aggregation + RDR.

All bulk assays (WGS/WGS-lr/WES) share ONE set of fixed bins tiled from ``segment.bed``
(the window BED). ``build_adaptive_bins`` closes a bb only when every tumor observation
meets ``min_snp_reads``, clustered by ``seg_id`` (the segment). Fixed-bin depth is
aggregated per assay onto the same bbs. Allele counts are aggregated per bb across all
samples; RDR is computed per assay, normalizing each tumor by the RDR base observation
named in its ``RDR_BASE_REP_ID`` (median-normalized when unset).

The allele matrices come as one joint set from phase_and_concat_bulk (read directly, no
union); depth/fixed-bin inputs stay per-assay (index-aligned to ``params.assay_types``).
Outputs go under ``bb_dir/MSR{msr}/bulk/``; matrix observations are the bulk samples.
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, log_hist, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd

from aggregation_utils import (
    build_adaptive_bins,
    stamp_bin_label,
    log_off_range_depth,
)
from range_utils import assign_pos_to_range
from feature_utils import (
    explode_feature_ids,
    stamp_bb_feature_ids,
    stamp_gene_clusters,
)
from io_utils import read_snp_mats, read_window_bed, write_bb_file
from matrix_utils import sum_features_to_bbs
from combine_counts_utils import (
    aggregate_bin_depth_to_bbs,
    build_rdr_base_map,
    compute_bb_rdr,
)
from phasing_utils import (
    detect_phase_flips,
    estimate_switchprobs_PS,
    estimate_switchprobs_cM,
    interp_cM_between_bbs,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_combine_counts import plot_rdr_baf, plot_rdr_baf_2d, plot_segmentation_qc

##################################################

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
dp_corrected_files = list(snakemake_handle.input["dp_corrected"])
bin_df_files = list(snakemake_handle.input["window_df"])
sample_file = snakemake_handle.input["sample_file"]
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])
genome_size = snakemake_handle.input["genome_size"]

# parameters
sample_id = snakemake_handle.params["sample_id"]
assay_types = list(snakemake_handle.params["assay_types"])

phase_flip_test = bool(snakemake_handle.params["phase_flip_test"])
phase_flip_epsilon = float(snakemake_handle.params["phase_flip_epsilon"])
phase_flip_alpha = float(snakemake_handle.params["phase_flip_alpha"])

gene_aware_binning = bool(snakemake_handle.params["gene_aware_binning"])
max_blocksize = int(snakemake_handle.params["max_blocksize"])
msr_list = [int(m) for m in snakemake_handle.params["min_snp_reads"]]
min_snp_per_bin = int(snakemake_handle.params["min_snp_per_bin"])
rdr_outlier_quantile = float(snakemake_handle.params["rdr_outlier_quantile"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])

# outputs
out_bb_file = list(snakemake_handle.output["bb_file"])
out_tot_mtx_bb = list(snakemake_handle.output["tot_mtx_bb"])
out_a_mtx_bb = list(snakemake_handle.output["a_mtx_bb"])
out_b_mtx_bb = list(snakemake_handle.output["b_mtx_bb"])
out_dp_mtx_bb = list(snakemake_handle.output["dp_mtx_bb"])
out_rdr_mtx_bb = list(snakemake_handle.output["rdr_mtx_bb"])
out_sample_file = list(snakemake_handle.output["sample_file"])
out_qc_pdf = list(snakemake_handle.output["qc_pdf"])

##################################################
# load inputs
sample_df = pd.read_table(sample_file)
bin_df, dp_bin_dfs = read_window_bed(bin_df_files)
snps, tot_mtx, a_mtx, b_mtx = read_snp_mats(
    snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp, mat_dtype=np.int32
)

num_datasets = len(sample_df)
dataset_assays = sample_df["assay_type"].to_numpy()
dataset_ids = sample_df["REP_ID"].tolist()
sample_types = sample_df["sample_type"].to_numpy()
assay2dataset_inds = {at: dataset_assays == at for at in assay_types}
empty = [at for at, inds in assay2dataset_inds.items() if not inds.any()]
assert not empty, f"joint sample sheet, no dataset for assay(s) {empty}"
tumor_obs = np.flatnonzero(sample_types == "tumor").tolist()
base_map = build_rdr_base_map(sample_df)
logging.info(
    f"combine_counts\n"
    f"sample_id={sample_id}\n"
    f"assay_types={assay_types}\n"
    f"#SNPs={len(snps)}\n"
    f"#windows={len(bin_df)}\n"
    f"#datasets={num_datasets}\n"
    f"#tumor_datasets={len(tumor_obs)}"
)

if "PS" not in snps.columns:
    snps["PS"] = 1
assert snps["PS"].notna().all(), "SNP file, `PS` column has NaNs"
cluster_cols = ["region_id", "seg_id", "PS"]

if phase_flip_test:
    snps["phase_cluster"] = detect_phase_flips(
        snps,
        a_mtx[:, tumor_obs],
        b_mtx[:, tumor_obs],
        cluster_cols=cluster_cols,
        tumor_sidx=0,
        epsilon=phase_flip_epsilon,
        alpha=phase_flip_alpha,
    )
    cluster_cols.append("phase_cluster")

tot_tumor = np.ascontiguousarray(tot_mtx[:, tumor_obs], dtype=np.float64)
# assigned once here; every sweep point below reuses it
snps["_orig_df_idx"] = np.arange(len(snps))
snps_binned, off_idx = assign_pos_to_range(snps, bin_df, ref_id="bin_id", dropna=True)
log_off_range_depth(off_idx, tot_tumor)
log_hist(
    snps_binned.groupby("bin_id").size().reindex(range(len(bin_df)), fill_value=0),
    "SNPs per fixed bin",
)
stamp_bin_label(bin_df, snps_binned, "PS", default=1)

if phase_flip_test:
    stamp_bin_label(bin_df, snps_binned, "phase_cluster", default=0)

if gene_aware_binning:
    stamp_gene_clusters(bin_df, snps_binned)

sample_labels = [
    f"{dataset_ids[i]} {dataset_assays[i]} {sample_types[i][0].upper()}"
    for i in range(num_datasets)
]
tumor_labels = [sample_labels[c] for c in tumor_obs]
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None

dp_corrected_list = [np.load(f)["mat"] for f in dp_corrected_files]

for msr, out_bb, out_tot, out_a, out_b, out_dp, out_rdr, out_samp, out_pdf in zip(
    msr_list,
    out_bb_file,
    out_tot_mtx_bb,
    out_a_mtx_bb,
    out_b_mtx_bb,
    out_dp_mtx_bb,
    out_rdr_mtx_bb,
    out_sample_file,
    out_qc_pdf,
):
    logging.info(f"===== MSR={msr} =====")
    min_snp_reads_vec = np.full(len(tumor_obs), msr, dtype=np.float64)
    bbs, snps_bb = build_adaptive_bins(
        bin_df,
        snps_binned.copy(),
        tot_tumor,
        min_snp_reads_vec,
        min_snp_per_bin,
        cluster_cols=cluster_cols,
        tumor_sidx=0,
        max_blocksize=max_blocksize,
        gene_aware=gene_aware_binning,
    )
    num_bbs = len(bbs)

    bb_ids = snps_bb["bb_id"].to_numpy()
    snp_orig_df_idx = snps_bb["_orig_df_idx"].to_numpy()
    a_mtx_bb = sum_features_to_bbs(a_mtx[snp_orig_df_idx], bb_ids, num_bbs)
    b_mtx_bb = sum_features_to_bbs(b_mtx[snp_orig_df_idx], bb_ids, num_bbs)
    tot_mtx_bb = sum_features_to_bbs(tot_mtx[snp_orig_df_idx], bb_ids, num_bbs)

    baf_mtx_bb = np.divide(
        b_mtx_bb,
        tot_mtx_bb,
        where=tot_mtx_bb > 0,
        out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
    )

    logging.info("aggregating corrected fixed-bin depth into bbs (per assay)")
    bb_dp, bb_bases = aggregate_bin_depth_to_bbs(
        assay2dataset_inds,
        bin_df,
        dp_bin_dfs,
        dp_corrected_list,
        num_bbs,
        num_datasets,
    )

    logging.info(
        f"compute bb RDR, {len(base_map)}/{len(tumor_obs)} tumors with RDR base"
    )
    bb_rdr = compute_bb_rdr(
        assay2dataset_inds,
        dp_bin_dfs,
        dp_corrected_list,
        bb_dp,
        tumor_obs,
        base_map,
        rdr_outlier_quantile,
        dataset_ids,
    )

    rdr_ylim = (np.round(np.nanquantile(bb_rdr, 0.99)).astype(int) + 1) * 1.1

    if gene_aware_binning:
        genic = explode_feature_ids(snps_bb, cols=["bb_id"])
        bb_gene_count = (
            genic.groupby("bb_id")["feature_id"]
            .nunique()
            .reindex(range(num_bbs))
            .fillna(0)
            .to_numpy()
        )
    else:
        bb_gene_count = None

    depth_tumor = bb_dp[:, tumor_obs]
    depth_normal = np.full_like(depth_tumor, np.nan, dtype=float)
    rdr_titles, rdr_norm_labels = [], []
    for j, c in enumerate(tumor_obs):
        title = f"{sample_id} - {sample_labels[c]}"
        if c in base_map:
            depth_normal[:, j] = bb_dp[:, base_map[c]]
            rdr_norm_labels.append("normal")
            rdr_titles.append(f"{title} / {sample_labels[base_map[c]]}")
        else:
            rdr_norm_labels.append("median")
            rdr_titles.append(title)

    baf_tumor = baf_mtx_bb[:, tumor_obs]

    with PdfPages(out_pdf) as pdf:
        plot_segmentation_qc(
            bbs,
            sample_df,
            bb_bases,
            b_mtx_bb,
            tot_mtx_bb,
            pdf=pdf,
            gene_count=bb_gene_count,
            sample_id=sample_id,
        )
        plot_rdr_baf(
            bbs,
            bb_rdr,
            baf_tumor,
            depth_tumor,
            depth_normal,
            rdr_titles,
            rdr_norm_labels,
            genome_size,
            out_pdf,
            feature_label="bb",
            rdr_ylim=rdr_ylim,
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
            pdf=pdf,
        )
        plot_rdr_baf_2d(
            bb_rdr,
            baf_mtx_bb[:, tumor_obs],
            tumor_labels,
            rdr_ylim=rdr_ylim,
            pdf=pdf,
        )
    logging.info(f"saved QC PDF to {out_pdf}")

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
        bb_dp = bb_dp[valid]
        bb_rdr = bb_rdr[valid]

    kept_bb_ids = np.where(~nan_mask)[0]
    old_to_new = {old: new for new, old in enumerate(kept_bb_ids)}
    snps_valid = snps_bb[snps_bb["bb_id"].isin(old_to_new)].copy()
    snps_valid["bb_id"] = snps_valid["bb_id"].map(old_to_new)
    bbs["bb_id"] = np.arange(len(bbs))

    logging.info("estimate bin-level switchprobs")
    if genetic_map is not None:
        dist_cms = interp_cM_between_bbs(
            bbs, snps_valid, genetic_map, bb_id_col="bb_id"
        )
        bbs["switchprobs"] = estimate_switchprobs_cM(
            dist_cms,
            nu=nu,
            min_switchprob=min_switchprob,
        )
    else:
        bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

    stamp_bb_feature_ids(bbs, snps_valid)
    write_bb_file(bbs, out_bb)
    np.savez_compressed(out_tot, mat=tot_mtx_bb)
    np.savez_compressed(out_a, mat=a_mtx_bb)
    np.savez_compressed(out_b, mat=b_mtx_bb)
    np.savez_compressed(out_dp, mat=bb_dp)
    np.savez_compressed(out_rdr, mat=bb_rdr)
    sample_df.to_csv(out_samp, sep="\t", index=False)
    logging.info(f"MSR={msr}: wrote {len(bbs)} bins")

# TODO: recommend a default MSR (e.g. elbow of lag-1 RDR/BAF dispersion vs #bins;
# see docs/combine_counts_pseudocode.md section 4) and record the pick.
logging.info("finished combine_counts.")
