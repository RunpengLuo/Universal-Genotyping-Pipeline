"""SNP-informed adaptive binning across all bulk assays + depth aggregation + RDR.

All bulk assays (WGS/WGS-lr/WES) share ONE set of fixed bins tiled from ``segment.bed``
(the window BED). ``build_adaptive_bins`` closes a bb only when every tumor observation
meets ``min_snp_reads``, clustered by ``seg_id`` (breakpoint chunk). Fixed-bin depth is
aggregated per assay onto the same bbs. Allele counts are aggregated per bb across all
samples; RDR is computed per assay, normalizing each tumor by the RDR base observation
named in its ``RDR_BASE_REP_ID`` (median-normalized when unset).

The allele matrices come as one joint set from phase_and_concat_bulk (read directly, no
union); depth/fixed-bin inputs stay per-assay (index-aligned to ``params.bulk_assays``).
Outputs go under ``bb_dir/MSR{msr}/bulk/``; matrix observations are the bulk samples.
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, maybe_path, sort_df_chr

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd

from aggregation_utils import (
    build_adaptive_bins,
    gene_cluster_labels,
    merge_feature_ids,
    assign_snps_to_bins,
)
from io_utils import read_snp_mats_bulk
from matrix_utils import sum_features_to_bbs
from combine_counts_utils import (
    aggregate_bin_depth_to_bbs,
    build_assay_obs_clusters,
    build_rdr_base_map,
    compute_bb_rdr,
)
from phasing_utils import (
    detect_phase_flips,
    estimate_switchprobs_PS,
    estimate_switchprobs_cM,
    interp_cM_between_bbs,
    setup_phase_clusters,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_combine_counts import plot_rdr_baf, plot_rdr_baf_2d, plot_segmentation_qc
from plot_utils import observation_order


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
blacklist_bed = maybe_path(snakemake_handle.input.get("blacklist_bed", None))
genome_size = snakemake_handle.input["genome_size"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
run_id = snakemake_handle.params["run_id"]
bulk_assays = list(snakemake_handle.params["bulk_assays"])
phase_flip_test = bool(snakemake_handle.params["phase_flip_test"])
phase_flip_epsilon = float(snakemake_handle.params["phase_flip_epsilon"])
phase_flip_alpha = float(snakemake_handle.params["phase_flip_alpha"])
gene_aware_binning_param = bool(snakemake_handle.params["gene_aware_binning"])
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

sample_df = pd.read_table(sample_file)
snps, tot_mtx, a_mtx, b_mtx = read_snp_mats_bulk(
    snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp
)
dp_corrected_list = [np.load(f)["mat"] for f in dp_corrected_files]
bin_df_list = [pd.read_table(f, sep="\t") for f in bin_df_files]
n_snps = len(snps)

sample_id = sample_df["SAMPLE_NAME"].iloc[0]
total_samples = len(sample_df)
logging.info(
    f"combine_counts: sample_id={sample_id}, bulk_assays={bulk_assays}; "
    f"{n_snps} SNPs x {total_samples} samples"
)

obs_assay = sample_df["assay_type"].tolist()
obs_repid = sample_df["REP_ID"].tolist()
assay_obs_clusters, tumor_obs_all = build_assay_obs_clusters(sample_df, bulk_assays)
base_map = build_rdr_base_map(sample_df)
logging.info(f"{total_samples} bulk samples, {len(tumor_obs_all)} tumor columns")

has_feature = "feature_id" in snps.columns
cluster_cols = setup_phase_clusters(snps)

if phase_flip_test:
    snps["phase_cluster"] = detect_phase_flips(
        snps,
        a_mtx[:, tumor_obs_all],
        b_mtx[:, tumor_obs_all],
        cluster_cols=cluster_cols,
        tumor_sidx=0,
        epsilon=phase_flip_epsilon,
        alpha=phase_flip_alpha,
    )
    cluster_cols.append("phase_cluster")

# one shared fixed-bin set: every bulk assay (WGS/WGS-lr/WES) tiles the same segment.bed
logging.info(f"fixed bins shared across assays {bulk_assays}")
_bcols = ["#CHR", "START", "END", "region_id"]
if all("seg_id" in w.columns for w in bin_df_list):
    _bcols.append("seg_id")
bin_df = pd.concat(
    [w[_bcols] for w in bin_df_list],
    ignore_index=True,
).drop_duplicates(["#CHR", "START", "END"])
bin_df = sort_df_chr(bin_df, ch="#CHR", pos="START").reset_index(drop=True)
if "seg_id" not in bin_df.columns:
    # no global BED seg_id on the fixed bins -> one seg per arm (== region_id partition)
    bin_df["seg_id"] = bin_df["region_id"]

gene_aware_binning = gene_aware_binning_param and has_feature
bin_df["bin_id"] = np.arange(len(bin_df))
tot_tumor = np.ascontiguousarray(tot_mtx[:, tumor_obs_all], dtype=np.float64)
# assigned once here; every sweep point below reuses it
snps_binned = assign_snps_to_bins(snps, bin_df, tot_tumor)
snps_per_bin = snps_binned.groupby("bin_id").size()
logging.info(
    f"SNPs per fixed bin: {len(snps_per_bin)}/{len(bin_df)} bins have SNPs, "
    f"mean={snps_per_bin.mean():.1f}, median={snps_per_bin.median():.1f}"
)
bin_ps = snps_binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
bin_df["PS"] = bin_df["bin_id"].map(bin_ps)
if bin_df["PS"].isna().any():
    bin_df["PS"] = bin_df["PS"].ffill()

if phase_flip_test:
    bin_pc = snps_binned.groupby("bin_id")["phase_cluster"].agg(
        lambda x: x.mode().iloc[0]
    )
    bin_df["phase_cluster"] = bin_df["bin_id"].map(bin_pc)
    if bin_df["phase_cluster"].isna().any():
        bin_df["phase_cluster"] = bin_df["phase_cluster"].ffill()

if gene_aware_binning:
    # glue each gene span into one cluster so a bb never splits a gene;
    # explode the ;-joined multi-gene feature_id so each gene gets its own span
    genic = snps_binned[
        snps_binned["feature_id"].notna() & (snps_binned["feature_id"] != "intergenic")
    ].copy()
    genic["feature_id"] = genic["feature_id"].str.split(";")
    genic = genic.explode("feature_id")
    genic = genic[genic["feature_id"] != "intergenic"]
    rng = genic.groupby("feature_id")["bin_id"].agg(["min", "max"])
    bin_df["gene_cluster"] = gene_cluster_labels(
        len(bin_df), zip(rng["min"].to_numpy(), rng["max"].to_numpy())
    )
    logging.info(
        f"gene-aware binning: {len(rng)} genes over {len(bin_df)} fixed bins -> "
        f"{bin_df['gene_cluster'].nunique()} gene/intergenic clusters (bbs never split a gene)"
    )

sample_labels = [f"{obs_assay[i]}:{obs_repid[i]}" for i in range(total_samples)]
tumor_labels = [sample_labels[c] for c in tumor_obs_all]
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None

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
    logging.info(f"===== binning MSR={msr} =====")
    min_snp_reads_vec = np.full(len(tumor_obs_all), msr, dtype=np.float64)
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
        assay_obs_clusters,
        bin_df,
        bin_df_list,
        dp_corrected_list,
        num_bbs,
        total_samples,
    )

    logging.info(
        f"compute bb RDR, {len(base_map)}/{len(tumor_obs_all)} tumors with RDR base"
    )
    bb_rdr = compute_bb_rdr(
        assay_obs_clusters,
        bin_df_list,
        dp_corrected_list,
        bb_dp,
        tumor_obs_all,
        base_map,
        rdr_outlier_quantile,
        obs_repid,
    )

    rdr_ylim = (np.round(np.nanquantile(bb_rdr, 0.99)).astype(int) + 1) * 1.1

    if gene_aware_binning and "feature_id" in snps_bb.columns:
        _genic = snps_bb[
            snps_bb["feature_id"].notna() & (snps_bb["feature_id"] != "intergenic")
        ][["bb_id", "feature_id"]].copy()
        _genic["feature_id"] = _genic["feature_id"].str.split(";")
        _genic = _genic.explode("feature_id")
        _genic = _genic[_genic["feature_id"] != "intergenic"]
        bb_gene_count = (
            _genic.groupby("bb_id")["feature_id"]
            .nunique()
            .reindex(range(num_bbs))
            .fillna(0)
            .to_numpy()
        )
    else:
        bb_gene_count = None

    depth_tumor = bb_dp[:, tumor_obs_all]
    depth_normal = np.full_like(depth_tumor, np.nan, dtype=float)
    rdr_titles, rdr_norm_labels = [], []
    for j, c in enumerate(tumor_obs_all):
        title = f"{sample_id} ({obs_assay[c]}) {obs_repid[c]} (T)"
        if c in base_map:
            depth_normal[:, j] = bb_dp[:, base_map[c]]
            rdr_norm_labels.append("normal")
            rdr_titles.append(f"{title} & {obs_repid[base_map[c]]} (N)")
        else:
            rdr_norm_labels.append("median")
            rdr_titles.append(title)

    # page order: assay (WGS<WGS-lr<WES) then dataset_id (all tumors here)
    t_order = observation_order(
        [obs_assay[c] for c in tumor_obs_all],
        ["tumor"] * len(tumor_obs_all),
        [obs_repid[c] for c in tumor_obs_all],
    )
    baf_tumor = baf_mtx_bb[:, tumor_obs_all]

    with PdfPages(out_pdf) as pdf:
        plot_segmentation_qc(
            bbs,
            sample_df,
            bb_bases,
            b_mtx_bb,
            tot_mtx_bb,
            pdf=pdf,
            gene_count=bb_gene_count,
        )
        plot_rdr_baf(
            bbs,
            bb_rdr[:, t_order],
            baf_tumor[:, t_order],
            depth_tumor[:, t_order],
            depth_normal[:, t_order],
            [rdr_titles[i] for i in t_order],
            [rdr_norm_labels[i] for i in t_order],
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
            baf_mtx_bb[:, tumor_obs_all],
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

    kept_bb_ids = np.where(~nan_mask)[0] if n_nan_rows > 0 else np.arange(num_bbs)
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

    bb_cols = ["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]
    if "feature_id" in snps_valid.columns:
        bbs["feature_id"] = (
            bbs["bb_id"]
            .map(snps_valid.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
            .fillna("intergenic")
        )
        bb_cols.append("feature_id")
    bbs[bb_cols].to_csv(out_bb, sep="\t", header=True, index=False)
    np.savez_compressed(out_tot, mat=tot_mtx_bb)
    np.savez_compressed(out_a, mat=a_mtx_bb)
    np.savez_compressed(out_b, mat=b_mtx_bb)
    np.savez_compressed(out_dp, mat=bb_dp)
    np.savez_compressed(out_rdr, mat=bb_rdr)
    sample_df.to_csv(out_samp, sep="\t", index=False)
    logging.info(f"MSR={msr}: wrote {len(bbs)} bins")

# TODO: recommend a default MSR (e.g. elbow of lag-1 RDR/BAF dispersion vs #bins;
# see docs/combine_counts_pseudocode.md section 4) and record the pick.
logging.info("finished combine_counts (all MSR).")
