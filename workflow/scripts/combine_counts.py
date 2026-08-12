"""Bulk: SNP-informed adaptive binning over all bulk assays, plus depth and RDR.

Last update: 2026-08-11

Inputs:
- allele_dir/snps.tsv.gz: the shared SNP set, matrix rows
- allele_dir/snp.{T,A,B}allele.npz: joint allele counts, samples as columns
- allele_dir/sample_ids.tsv: one row per matrix column
- aux_dir/windows.bed.gz: the shared fixed bins, matrix rows of the depth
- pileup_dir/bulk/window.dp.npz: bias-corrected depth, windows x bulk datasets
- aux_dir/segment.bed: region_id and seg_id cluster keys
- phase_dir/genetic_map.tsv.gz: optional, for cM-based switch probabilities
- blacklist_bed, genome_size: QC plot shading and axis
Outputs:
- bb_dir/MSR{msr}/bulk/bb.tsv.gz: bb definitions, one row each
- bb_dir/MSR{msr}/bulk/bb.{T,A,B}allele.npz: per-bb phased allele counts
- bb_dir/MSR{msr}/bulk/bb.depth.npz: per-bb mean depth per dataset
- bb_dir/MSR{msr}/bulk/bb.rdr.npz: per-bb RDR per tumor
- bb_dir/MSR{msr}/bulk/sample_ids.tsv: one row per matrix column
- qc_dir/combine_counts.bulk.MSR{msr}.pdf: segmentation, RDR/BAF and 2D QC
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, log_hist, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd

from segmentation_utils import build_adaptive_bins, sum_features_to_bbs
from range_utils import assign_pos_to_range, merge_ranges_to_clusters
from feature_utils import explode_feature_ids, merge_feature_ids
from io_utils import read_snp_mats, read_window_bed, write_bb_file
from combine_counts_utils import summarize_read_depth_bb, summarize_rdr_bb
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
dp_corrected_file = snakemake_handle.input["dp_corrected"]
window_bed = snakemake_handle.input["window_bed"]
sample_file = snakemake_handle.input["sample_file"]
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])
genome_size = snakemake_handle.input["genome_size"]

# parameters
sample_id = snakemake_handle.params["sample_id"]
assay_types = list(snakemake_handle.params["assay_types"])
dp_dataset_ids = list(snakemake_handle.params["dataset_ids"])
dp_dataset_assays = list(snakemake_handle.params["dataset_assays"])
chroms = list(snakemake_handle.params["chroms"])

phase_flip_test = bool(snakemake_handle.params["phase_flip_test"])
phase_flip_epsilon = float(snakemake_handle.params["phase_flip_epsilon"])
phase_flip_alpha = float(snakemake_handle.params["phase_flip_alpha"])

gene_aware_binning = bool(snakemake_handle.params["gene_aware_binning"])
max_blocksize = int(snakemake_handle.params["max_blocksize"])
msr_list = [int(m) for m in snakemake_handle.params["min_snp_reads"]]
min_snp_per_bin = int(snakemake_handle.params["min_snp_per_bin"])
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
bin_df = read_window_bed(window_bed, chroms=chroms)
snps, tot_mtx, a_mtx, b_mtx = read_snp_mats(
    snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp, mat_dtype=np.int32
)
dp_corrected = np.load(dp_corrected_file)["mat"]
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None

##################################################
# observation order, labels and RDR bases
num_datasets = len(sample_df)
dataset_assays = sample_df["assay_type"].to_numpy()
dataset_ids = sample_df["dataset_id"].tolist()
sample_types = sample_df["sample_type"].to_numpy()
tumor_dataset_indices = np.flatnonzero(sample_types == "tumor").tolist()
tumor_dataset_ids = [dataset_ids[c] for c in tumor_dataset_indices]
tumor_assays = [dataset_assays[c] for c in tumor_dataset_indices]

assert dp_dataset_ids == dataset_ids and dp_dataset_assays == list(dataset_assays), (
    f"depth columns {list(zip(dp_dataset_ids, dp_dataset_assays))} do not match "
    f"{sample_file}: {list(zip(dataset_ids, dataset_assays))}"
)
assert dp_corrected.shape == (len(bin_df), num_datasets), (
    f"{dp_corrected_file}: shape {dp_corrected.shape}, expected "
    f"({len(bin_df)}, {num_datasets}) from {window_bed}"
)

get_rdr_base_dataset_id = {}
if "rdr_base_dataset_id" in sample_df.columns:
    obs_of = {rid: i for i, rid in enumerate(dataset_ids)}
    for i in tumor_dataset_indices:
        base = sample_df["rdr_base_dataset_id"].iloc[i]
        if pd.notna(base):
            get_rdr_base_dataset_id[i] = obs_of[base]
tumor_base_dataset_ids = [
    dataset_ids[get_rdr_base_dataset_id[c]] if c in get_rdr_base_dataset_id else None
    for c in tumor_dataset_indices
]

logging.info(
    f"combine_counts\n"
    f"sample_id={sample_id}\n"
    f"assay_types={assay_types}\n"
    f"#SNPs={len(snps)}\n"
    f"#windows={len(bin_df)}\n"
    f"#datasets={num_datasets}\n"
    f"#tumor_datasets={len(tumor_dataset_indices)}"
)

##################################################
# assign SNPs to windows, then drop the misses from the matrices too
tot_tumor = np.ascontiguousarray(tot_mtx[:, tumor_dataset_indices], dtype=np.float64)
snps_binned, off_idx = assign_pos_to_range(snps, bin_df, ref_id="bin_id", dropna=True)
if len(off_idx):
    log_hist(tot_tumor[off_idx].sum(axis=1), "depth of SNPs outside every bin")
log_hist(
    snps_binned.groupby("bin_id").size().reindex(range(len(bin_df)), fill_value=0),
    "SNPs per window",
)
keep_snps = np.ones(len(snps), dtype=bool)
keep_snps[off_idx] = False
tot_mtx, a_mtx, b_mtx = tot_mtx[keep_snps], a_mtx[keep_snps], b_mtx[keep_snps]
tot_tumor = np.ascontiguousarray(tot_tumor[keep_snps])

##################################################
# adaptive segmentation bounderies
cluster_cols = ["region_id", "seg_id"]

if "PS" in snps.columns:
    assert snps["PS"].notna().all(), "SNP file, `PS` column has NaNs"
    cluster_cols.append("PS")
    modal = snps_binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
    bin_df["PS"] = bin_df["bin_id"].map(modal).ffill().bfill().fillna(1)

if phase_flip_test:
    snps_binned["phase_cluster"] = detect_phase_flips(
        snps_binned,
        a_mtx[:, tumor_dataset_indices],
        b_mtx[:, tumor_dataset_indices],
        cluster_cols=cluster_cols,
        epsilon=phase_flip_epsilon,
        alpha=phase_flip_alpha,
    )
    cluster_cols.append("phase_cluster")
    modal = snps_binned.groupby("bin_id")["phase_cluster"].agg(
        lambda x: x.mode().iloc[0]
    )
    bin_df["phase_cluster"] = bin_df["bin_id"].map(modal).ffill().bfill().fillna(0)

if gene_aware_binning:
    gene_spans = (
        explode_feature_ids(snps_binned, cols=["bin_id"])
        .groupby("feature_id")["bin_id"]
        .agg(["min", "max"])
    )
    bin_df["gene_cluster"] = merge_ranges_to_clusters(
        len(bin_df), zip(gene_spans["min"].to_numpy(), gene_spans["max"].to_numpy() + 1)
    )

##################################################
# one adaptive binning per min_snp_reads, on the shared fixed bins
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
    min_snp_reads_vec = np.full(len(tumor_dataset_indices), msr, dtype=np.float64)
    bbs, snps_bb = build_adaptive_bins(
        bin_df,
        snps_binned,
        tot_tumor,
        min_snp_reads_vec,
        min_snp_per_bin,
        cluster_cols=cluster_cols,
        max_blocksize=max_blocksize,
        gene_aware=gene_aware_binning,
    )
    num_bbs = len(bbs)

    bb_ids = snps_bb["bb_id"].to_numpy()
    a_mtx_bb = sum_features_to_bbs(a_mtx, bb_ids, num_bbs)
    b_mtx_bb = sum_features_to_bbs(b_mtx, bb_ids, num_bbs)
    tot_mtx_bb = sum_features_to_bbs(tot_mtx, bb_ids, num_bbs)

    baf_mtx_bb = np.divide(
        b_mtx_bb,
        tot_mtx_bb,
        where=tot_mtx_bb > 0,
        out=np.full_like(b_mtx_bb, np.nan, dtype=np.float32),
    )

    logging.info("aggregating corrected fixed-bin depth into bbs")
    bb_dp, bb_bases = summarize_read_depth_bb(bin_df, dp_corrected, num_bbs)

    logging.info(
        f"compute bb RDR, {len(get_rdr_base_dataset_id)}/{len(tumor_dataset_indices)} tumors with RDR base"
    )
    bb_rdr = summarize_rdr_bb(
        bin_df,
        dp_corrected,
        bb_dp,
        tumor_dataset_indices,
        get_rdr_base_dataset_id,
        dataset_ids,
    )

    # drop bbs carrying a NaN, before the QC plots, so the PDF shows what is written
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
        bb_bases = bb_bases[valid]
        bb_rdr = bb_rdr[valid]
    num_bbs = len(bbs)

    kept_bb_ids = np.where(~nan_mask)[0]
    old_to_new = {old: new for new, old in enumerate(kept_bb_ids)}
    snps_valid = snps_bb[snps_bb["bb_id"].isin(old_to_new)].copy()
    snps_valid["bb_id"] = snps_valid["bb_id"].map(old_to_new)
    bbs["bb_id"] = np.arange(num_bbs)

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

    bbs["feature_id"] = (
        bbs["bb_id"]
        .map(snps_valid.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )

    depth_tumor = bb_dp[:, tumor_dataset_indices]
    depth_normal = np.full_like(depth_tumor, np.nan, dtype=float)
    for j, c in enumerate(tumor_dataset_indices):
        if c in get_rdr_base_dataset_id:
            depth_normal[:, j] = bb_dp[:, get_rdr_base_dataset_id[c]]

    baf_tumor = baf_mtx_bb[:, tumor_dataset_indices]

    with PdfPages(out_pdf) as pdf:
        plot_segmentation_qc(
            bbs,
            sample_df,
            bb_bases,
            b_mtx_bb,
            tot_mtx_bb,
            pdf=pdf,
            sample_id=sample_id,
        )
        plot_rdr_baf(
            bbs,
            bb_rdr,
            baf_tumor,
            depth_tumor,
            depth_normal,
            sample_id,
            tumor_dataset_ids,
            tumor_assays,
            tumor_base_dataset_ids,
            genome_size,
            out_pdf,
            feature_label="bb",
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
            pdf=pdf,
        )
        plot_rdr_baf_2d(
            bb_rdr,
            baf_tumor,
            sample_id,
            tumor_dataset_ids,
            tumor_assays,
            pdf=pdf,
        )
    logging.info(f"saved QC PDF to {out_pdf}")

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
