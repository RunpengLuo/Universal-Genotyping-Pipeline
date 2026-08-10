"""Joint per-sample adaptive binning for non-bulk assays.

All non-bulk assays present in the sample (e.g. scRNA + scATAC for multiome) are segmented
onto ONE shared set of bbs. Each (replicate x assay) is pseudobulked into one observation;
``build_adaptive_bins`` then requires ``min_snp_reads`` in EVERY observation, so bbs jointly
satisfy every (dataset_id, assay) -- exactly the bulk multi-sample pattern (see combine_counts.py),
but with single cells pseudobulked per replicate first.

The fixed bins are the window BED, the same grid bulk bins on: windows exist where no SNP
does, so a SNP-free segment still yields bbs carrying Xcount, and no window lies in a
blacklist hole. The two count types reach a bb differently, on purpose: scATAC fragments
are routed through the windows (a fragment in a hole hits no window and is dropped), while
scRNA/VISIUM genes go to the bb hull, since a gene is indivisible and would land in an
arbitrary window otherwise.

All outputs live under ``bb_dir/MSR{msr}/{assay}/``: the shared ``bb.tsv.gz`` and combined
``sample_ids.tsv`` (duplicated per assay) plus the per-assay (n_bins x n_cells) matrices
(``bb.{T,A,B}allele.npz``, ``multi_snp.*``, ``barcodes*``). Input for HATCHet3 and CalicoST.
"""

import logging
import shutil


snakemake_handle = snakemake

from utils import log_hist, maybe_path, set_omp_threads, setup_logging

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from scipy.sparse import save_npz

from const import ASSAY_TYPE2MODALITY
from io_utils import (
    read_barcodes,
    read_full_barcodes,
    read_snp_mats,
    read_window_bed,
    write_bb_file,
)
from combine_counts_utils import (
    build_union_snps,
    observation_cluster_ids,
)
from phasing_utils import (
    estimate_switchprobs_PS,
    estimate_switchprobs_cM,
    interp_cM_between_bbs,
)
from matrix_utils import sum_features_to_bbs, sum_observations_to_pseudobulk
from feature_utils import (
    explode_feature_ids,
    merge_feature_ids,
    sum_atac_fragments_to_bins,
    sum_umis_to_bins,
)
from aggregation_utils import build_adaptive_bins
from range_utils import assign_pos_to_range, merge_ranges_to_clusters
from plot_alleles import plot_allele_freqs
from matplotlib.backends.backend_pdf import PdfPages

##################################################

# inputs
snp_info_files = list(snakemake_handle.input["snp_info"])
tot_mtx_snp_files = list(snakemake_handle.input["tot_mtx_snp"])
a_mtx_snp_files = list(snakemake_handle.input["a_mtx_snp"])
b_mtx_snp_files = list(snakemake_handle.input["b_mtx_snp"])
sample_files = list(snakemake_handle.input["sample_file"])
barcode_files = list(snakemake_handle.input["all_barcodes"])
barcode_full_files = list(snakemake_handle.input["barcodes_full"])
frag_files = list(snakemake_handle.input["frag_files"])
h5ad_files = list(snakemake_handle.input["h5ad_files"])
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
window_bed = snakemake_handle.input["window_bed"]
genome_size = snakemake_handle.input["genome_size"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
run_id = snakemake_handle.params["run_id"]
assay_types = list(snakemake_handle.params["assay_types"])
chroms = list(snakemake_handle.params["chroms"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])
nsnp_multi = int(snakemake_handle.params["nsnp_multi"])
msr_list = [int(m) for m in snakemake_handle.params["min_snp_reads"]]
min_snp_per_bin = int(snakemake_handle.params["min_snp_per_bin"])
gene_aware_binning = bool(snakemake_handle.params["gene_aware_binning"])

# outputs
out_bb_file = list(snakemake_handle.output["bb_file"])
out_sample_file = list(snakemake_handle.output["sample_file"])
out_tot_mtx_bb = list(snakemake_handle.output["tot_mtx_bb"])
out_a_mtx_bb = list(snakemake_handle.output["a_mtx_bb"])
out_b_mtx_bb = list(snakemake_handle.output["b_mtx_bb"])
out_multi_snp_file = list(snakemake_handle.output["multi_snp_file"])
out_tot_mtx_multi = list(snakemake_handle.output["tot_mtx_multi"])
out_a_mtx_multi = list(snakemake_handle.output["a_mtx_multi"])
out_b_mtx_multi = list(snakemake_handle.output["b_mtx_multi"])
out_all_barcodes = list(snakemake_handle.output["all_barcodes"])
out_barcodes_full = list(snakemake_handle.output["barcodes_full"])
out_x_count = list(snakemake_handle.output["x_count"])
out_qc_pdf = list(snakemake_handle.output["qc_pdf"])

n_assays = len(assay_types)

##################################################
# load per-assay inputs
sample_ids_list = [pd.read_table(f) for f in sample_files]
snps_list, tot_mtx_snp_list, a_mtx_snp_list, b_mtx_snp_list = (
    list(mats)
    for mats in zip(
        *[
            read_snp_mats(*files)
            for files in zip(
                snp_info_files, tot_mtx_snp_files, a_mtx_snp_files, b_mtx_snp_files
            )
        ]
    )
)
dataset_ids_list = [s["REP_ID"].tolist() for s in sample_ids_list]
sample_labels_list = [
    [
        f"{dataset_id} {assay_types[k]} {str(sample_type)[0].upper()}"
        for dataset_id, sample_type in zip(
            dataset_ids_list[k], sample_ids_list[k]["sample_type"]
        )
    ]
    for k in range(n_assays)
]
cell_dataset_idx_list = [
    observation_cluster_ids(read_full_barcodes(bc_full), dataset_ids_list[k])
    for k, bc_full in enumerate(barcode_full_files)
]

logging.info(
    f"combine_counts_nonbulk\n"
    f"sample_id={sample_id}\n"
    f"assay_types={assay_types}\n"
    f"#SNPs(per assay)={[len(s) for s in snps_list]}\n"
    f"#datasets(per assay)={[len(d) for d in dataset_ids_list]}"
)

##################################################
# 1. shared SNP set (union across assays)
snps = build_union_snps(snps_list)
n_snps = len(snps)
logging.info(f"shared SNP set (union): {n_snps} SNPs across {n_assays} assays")

if "PS" not in snps.columns:
    snps["PS"] = 1
assert snps["PS"].notna().all(), "SNP file, `PS` column has NaNs"
cluster_cols = ["region_id", "seg_id", "PS"]
logging.info(f"gene_aware_binning={gene_aware_binning}")

##################################################
# 2. per-(replicate x assay) pseudobulk scattered onto the shared SNP set
total_cols = sum(len(r) for r in dataset_ids_list)
tot_pb = np.zeros((n_snps, total_cols), dtype=np.float64)
tot_pb_list = []  # per-assay (n_snps_k x n_datasets_k) pseudobulk, in that assay's SNP order
offset = 0
for k in range(n_assays):
    n_datasets_k = len(dataset_ids_list[k])
    tot_pb_k = sum_observations_to_pseudobulk(
        tot_mtx_snp_list[k], cell_dataset_idx_list[k], n_datasets_k
    )
    tot_pb_list.append(tot_pb_k)
    shared_snp_ids = (
        snps_list[k][["#CHR", "POS0"]]
        .merge(snps[["#CHR", "POS0", "snp_id"]], on=["#CHR", "POS0"], how="left")[
            "snp_id"
        ]
        .to_numpy()
        .astype(np.int64)
    )
    tot_pb[shared_snp_ids, offset : offset + tot_pb_k.shape[1]] = tot_pb_k
    offset += n_datasets_k
logging.info(f"binning on {total_cols} (replicate x assay) pseudobulk observations")

##################################################
# 3. shared precompute (genetic map, sample sheet, fixed bins, multi-SNP clusters)
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None
rna_assay_types = [at for at in assay_types if ASSAY_TYPE2MODALITY[at] == "RNA"]
assert len(h5ad_files) == len(rna_assay_types), (
    f"h5ad_files, {len(h5ad_files)} files for {len(rna_assay_types)} RNA assays"
)
h5ad_by_assay = dict(zip(rna_assay_types, h5ad_files))
# combined sample sheet: one entry per (replicate x assay) observation
joint_sids = pd.concat(
    [sample_ids_list[k].assign(assay_type=assay_types[k]) for k in range(n_assays)],
    ignore_index=True,
)
# fixed bins: the window BED, the same grid bulk bins on. Tiled per segment row, so a
# window never spans two segments and none lies in a blacklist hole. Windows exist where
# no SNP does, so a SNP-free segment still yields bbs carrying Xcount.
bin_df, _ = read_window_bed(window_bed, chroms=chroms)
logging.info(f"fixed bins: {len(bin_df)} windows from {window_bed}")

tot_pb_cont = np.ascontiguousarray(tot_pb)
# assigned once here; every MSR below reuses it
snps_binned, off_idx = assign_pos_to_range(snps, bin_df, ref_id="bin_id", dropna=True)
if len(off_idx):
    log_hist(tot_pb_cont[off_idx].sum(axis=1), "depth of SNPs outside every bin")
keep_snps = np.ones(len(snps), dtype=bool)
keep_snps[off_idx] = False
tot_pb_cont = np.ascontiguousarray(tot_pb_cont[keep_snps])
modal = snps_binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
bin_df["PS"] = bin_df["bin_id"].map(modal).ffill().bfill().fillna(1)
if gene_aware_binning:
    gene_spans = (
        explode_feature_ids(snps_binned, cols=["bin_id"])
        .groupby("feature_id")["bin_id"]
        .agg(["min", "max"])
    )
    bin_df["gene_cluster"] = merge_ranges_to_clusters(
        len(bin_df), zip(gene_spans["min"].to_numpy(), gene_spans["max"].to_numpy() + 1)
    )
    logging.info(
        f"gene-aware binning: {len(gene_spans)} genes over {len(bin_df)} fixed bins -> "
        f"{bin_df['gene_cluster'].nunique()} clusters"
    )

# per-assay multi-SNP pre-grouping (diagnostic; every nsnp_multi SNPs)
multi_cache = []
for k in range(n_assays):
    snps_k = snps_list[k].copy()
    if "PS" not in snps_k.columns:
        snps_k["PS"] = 1
    k_cols = ["#CHR", "START", "END", "region_id"]
    if "seg_id" in snps_k.columns:
        k_cols.append("seg_id")
    bins_k = snps_k[k_cols].copy()
    bins_k["bin_id"] = np.arange(len(bins_k))
    snps_k, off_k = assign_pos_to_range(snps_k, bins_k, ref_id="bin_id", dropna=True)
    keep_k = np.ones(len(bins_k), dtype=bool)
    keep_k[off_k] = False
    if len(off_k):
        log_hist(tot_pb_list[k][off_k].sum(axis=1), "depth of SNPs outside every bin")
    multi_snps, snps_multi = build_adaptive_bins(
        bins_k,
        snps_k,
        np.ascontiguousarray(tot_pb_list[k][keep_k]),
        0,
        nsnp_multi,
        cluster_cols=[c for c in ("region_id", "seg_id") if c in bins_k.columns],
        max_blocksize=0,
        gene_aware=False,
    )
    multi_ids = snps_multi["bb_id"].to_numpy()
    n_multi = len(multi_snps)
    tot_multi = sum_features_to_bbs(tot_mtx_snp_list[k][keep_k], multi_ids, n_multi)
    a_multi = sum_features_to_bbs(a_mtx_snp_list[k][keep_k], multi_ids, n_multi)
    b_multi = sum_features_to_bbs(b_mtx_snp_list[k][keep_k], multi_ids, n_multi)
    if genetic_map is not None:
        dist_cms_multi = interp_cM_between_bbs(
            multi_snps, snps_multi, genetic_map, bb_id_col="bb_id"
        )
        multi_snps["switchprobs"] = estimate_switchprobs_cM(
            dist_cms_multi, nu=nu, min_switchprob=min_switchprob
        )
    else:
        multi_snps["switchprobs"] = estimate_switchprobs_PS(multi_snps, switchprob_ps)
    multi_cache.append(
        {
            "df": multi_snps.rename(columns={"bb_id": "multi_id"}),
            "tot": tot_multi,
            "a": a_multi,
            "b": b_multi,
        }
    )

##################################################
# 4. per-MSR joint segmentation + per-assay outputs
n_msr = len(msr_list)
for j, min_snp_reads in enumerate(msr_list):
    logging.info(f"===== MSR={min_snp_reads} =====")
    bbs, snps_bb = build_adaptive_bins(
        bin_df,
        snps_binned.copy(),
        tot_pb_cont,
        min_snp_reads,
        min_snp_per_bin,
        cluster_cols=cluster_cols,
        max_blocksize=0,
        gene_aware=gene_aware_binning,
    )
    num_bbs = len(bbs)
    assert bin_df["bb_id"].between(0, num_bbs - 1).all(), (
        "fixed bins, some did not land in a bb (null cluster key)"
    )
    n_empty = int((bbs["#SNPS"] == 0).sum())
    logging.info(
        f"{num_bbs} bbs from {len(bin_df)} windows; {n_empty} carry no SNP "
        "(Xcount only, all-zero allele rows)"
    )

    if genetic_map is not None:
        dist_cms = interp_cM_between_bbs(bbs, snps_bb, genetic_map, bb_id_col="bb_id")
        bbs["switchprobs"] = estimate_switchprobs_cM(
            dist_cms, nu=nu, min_switchprob=min_switchprob
        )
    else:
        bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

    bbs["feature_id"] = (
        bbs["bb_id"]
        .map(snps_bb.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )

    # union SNP -> shared bb_id map, plus the two frames counts are assigned through
    bb_of_snp = snps_bb[["#CHR", "POS0", "bb_id"]]
    # ATAC: the windows, stamped with their owning bb. A fragment in a blacklist hole or
    # between segments hits no window and is dropped, where the bb hull would swallow it.
    window_bb_ranges = bin_df[["#CHR", "START", "END", "bb_id"]].copy()
    # RNA: the bb hulls. A gene is never split, and against 1 kb windows every window
    # inside a gene ties on overlap, so gene -> window would pick one arbitrarily.
    bb_ranges = bbs[["#CHR", "START", "END", "bb_id"]].copy()

    for k in range(n_assays):
        assay = assay_types[k]
        # outputs are expanded assay-major over msr_list, so this is (k, j)
        idx = k * n_msr + j

        write_bb_file(bbs, out_bb_file[idx])
        joint_sids.to_csv(out_sample_file[idx], sep="\t", index=False)

        # map this assay's SNPs to the shared bbs
        snp_bbs = snps_list[k][["#CHR", "POS0"]].merge(
            bb_of_snp, on=["#CHR", "POS0"], how="left"
        )
        keep = snp_bbs["bb_id"].notna().to_numpy()
        snp_rows_k = np.where(keep)[0]
        bb_ids_k = snp_bbs.loc[keep, "bb_id"].to_numpy().astype(np.int64)

        tot_bb = sum_features_to_bbs(tot_mtx_snp_list[k][snp_rows_k], bb_ids_k, num_bbs)
        a_bb = sum_features_to_bbs(a_mtx_snp_list[k][snp_rows_k], bb_ids_k, num_bbs)
        b_bb = sum_features_to_bbs(b_mtx_snp_list[k][snp_rows_k], bb_ids_k, num_bbs)
        save_npz(out_tot_mtx_bb[idx], tot_bb)
        save_npz(out_a_mtx_bb[idx], a_bb)
        save_npz(out_b_mtx_bb[idx], b_bb)

        # per-cell Xcount per bb bin: scATAC from raw fragments, RNA from the h5ad
        if assay == "scATAC":
            assert len(frag_files) == len(dataset_ids_list[k]), (
                f"frag_files, {len(frag_files)} files for "
                f"{len(dataset_ids_list[k])} scATAC datasets"
            )
            x_count = sum_atac_fragments_to_bins(
                frag_files,
                dataset_ids_list[k],
                read_full_barcodes(barcode_full_files[k]),
                window_bb_ranges,
                num_bbs,
            )
            save_npz(out_x_count[idx], x_count)
            logging.info(
                f"{assay} MSR={min_snp_reads} Xcount (fragments): shape={x_count.shape}, nnz={x_count.nnz}"
            )
        elif assay in h5ad_by_assay:
            x_count = sum_umis_to_bins(
                h5ad_by_assay[assay],
                read_barcodes(barcode_files[k]),
                bb_ranges,
                num_bbs,
                assay,
            )
            save_npz(out_x_count[idx], x_count)
            logging.info(
                f"{assay} MSR={min_snp_reads} Xcount (h5ad): shape={x_count.shape}, nnz={x_count.nnz}"
            )

        multi = multi_cache[k]
        multi["df"].to_csv(out_multi_snp_file[idx], sep="\t", header=True, index=False)
        save_npz(out_tot_mtx_multi[idx], multi["tot"])
        save_npz(out_a_mtx_multi[idx], multi["a"])
        save_npz(out_b_mtx_multi[idx], multi["b"])

        with PdfPages(out_qc_pdf[idx]) as pdf:
            plot_allele_freqs(
                bbs,
                sample_labels_list[k],
                tot_bb,
                b_bb,
                genome_size,
                qc_dir,
                apply_pseudobulk=True,
                cell_dataset_ids=cell_dataset_idx_list[k],
                allele="B",
                feature_label="bb",
                run_id=f"{assay}.MSR{min_snp_reads}.{run_id}",
                name_prefix="combine_counts",
                sample_id=sample_id,
                pdf=pdf,
            )
            plot_allele_freqs(
                multi["df"],
                sample_labels_list[k],
                multi["tot"],
                multi["b"],
                genome_size,
                qc_dir,
                apply_pseudobulk=True,
                cell_dataset_ids=cell_dataset_idx_list[k],
                allele="B",
                feature_label="multi-snp",
                run_id=f"{assay}.MSR{min_snp_reads}.{run_id}",
                name_prefix="combine_counts",
                sample_id=sample_id,
                pdf=pdf,
            )

        shutil.copy2(barcode_files[k], out_all_barcodes[idx])
        shutil.copy2(barcode_full_files[k], out_barcodes_full[idx])

# TODO: recommend a default MSR (elbow of lag-1 dispersion vs #bins; see
# docs/combine_counts_pseudocode.md section 4) and record the pick.
logging.info("finished combine_counts_nonbulk.")
