"""Aggregate SNP-level and gene-level matrices onto pre-computed bbs.

Input for Copy-typing.
"""

import logging
import shutil


snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, sort_df_chr

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from scipy.sparse import issparse, load_npz, save_npz
import scanpy as sc
from const import ASSAY_TYPE2MODALITY, BULK_ASSAYS
from io_utils import read_barcodes, read_full_barcodes
from combine_counts_utils import observation_cluster_ids
from matrix_utils import sum_features_to_bbs
from range_utils import assign_pos_to_range
from feature_utils import (
    assign_features_to_ranges,
    merge_feature_ids,
    sum_atac_fragments_to_bins,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_alleles import plot_allele_freqs

COUNT_DTYPE = np.int32


def _sparsity(X):
    """Return fraction of zero entries."""
    size = X.shape[0] * X.shape[1]
    if size == 0:
        return 0.0
    nnz = X.nnz if issparse(X) else np.count_nonzero(X)
    return 1.0 - nnz / size


# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
h5ad_file = snakemake_handle.input["h5ad_file"]
frag_files = list(snakemake_handle.input["frag_files"])
all_barcodes = snakemake_handle.input["all_barcodes"]
barcodes_full_path = snakemake_handle.input["barcodes_full"]
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]
sample_file = snakemake_handle.input["sample_file"]
bb_file = snakemake_handle.input["bb_file"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
assay_type = snakemake_handle.params["assay_type"]
run_id = getattr(snakemake_handle.params, "run_id", "")

# outputs
out_x_count = snakemake_handle.output["x_count"]
out_tot_mtx_bb = snakemake_handle.output["tot_mtx_bb"]
out_a_mtx_bb = snakemake_handle.output["a_mtx_bb"]
out_b_mtx_bb = snakemake_handle.output["b_mtx_bb"]
cnv_segments = snakemake_handle.output["cnv_segments"]
barcodes_out = snakemake_handle.output["barcodes_out"]
barcodes_full_out = snakemake_handle.output["barcodes_full_out"]
out_sample_file = snakemake_handle.output["sample_file"]


sample_df = pd.read_table(sample_file)
dataset_ids = sample_df["REP_ID"].tolist()

is_rna_assay = ASSAY_TYPE2MODALITY[assay_type] == "RNA"
assert assay_type not in BULK_ASSAYS, (
    f"copytyping_preprocess, bulk assay unsupported: {assay_type}"
)

cell_rep_ids = observation_cluster_ids(
    read_full_barcodes(barcodes_full_path), dataset_ids
)

logging.info(f"cnv segmentation, sample_id={sample_id}, assay_type={assay_type}")
logging.info(f"dataset_ids={dataset_ids}")
snps = pd.read_table(snp_info, sep="\t")

tot_mtx = load_npz(tot_mtx_snp)
a_mtx = load_npz(a_mtx_snp)
b_mtx = load_npz(b_mtx_snp)

bb_df = pd.read_table(bb_file, sep="\t")
bb_df = sort_df_chr(bb_df, pos="START")
bb_df["bb_id"] = np.arange(len(bb_df))
num_bbs = len(bb_df)
logging.info(f"#bbs={num_bbs}")

snps["RAW_SNP_DF_IDX"] = np.arange(len(snps))
logging.info(f"#{assay_type}-SNP (raw)={len(snps)}")
snps, _ = assign_pos_to_range(snps, bb_df, ref_id="bb_id", dropna=True)
logging.info(f"#{assay_type}-SNP (remain)={len(snps)}")
bb_df["#SNPS"] = bb_df["bb_id"].map(snps["bb_id"].value_counts()).fillna(0).astype(int)

if "feature_id" in snps.columns:
    bb_df["feature_id"] = (
        bb_df["bb_id"]
        .map(snps.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )

raw_snp_df_idx = snps["RAW_SNP_DF_IDX"].to_numpy()
tot_mtx = tot_mtx[raw_snp_df_idx, :]
a_mtx = a_mtx[raw_snp_df_idx, :]
b_mtx = b_mtx[raw_snp_df_idx, :]

logging.info(
    f"SNP-level matrices: shape={tot_mtx.shape}, "
    f"tot sparsity={_sparsity(tot_mtx):.4f}, "
    f"A sparsity={_sparsity(a_mtx):.4f}, "
    f"B sparsity={_sparsity(b_mtx):.4f}"
)

bb_ids = snps["bb_id"].to_numpy()
tot_mtx_bb = sum_features_to_bbs(tot_mtx, bb_ids, num_bbs)
a_mtx_bb = sum_features_to_bbs(a_mtx, bb_ids, num_bbs)
b_mtx_bb = sum_features_to_bbs(b_mtx, bb_ids, num_bbs)
assert tot_mtx_bb.shape[0] == num_bbs

logging.info(
    f"bb-level matrices: shape={tot_mtx_bb.shape}, "
    f"T sparsity={_sparsity(tot_mtx_bb):.4f}, "
    f"A sparsity={_sparsity(a_mtx_bb):.4f}, "
    f"B sparsity={_sparsity(b_mtx_bb):.4f}"
)

pdf_path = snakemake_handle.output["qc_pdf"]
with PdfPages(pdf_path) as pdf:
    plot_allele_freqs(
        snps,
        dataset_ids,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_rep_ids=cell_rep_ids,
        allele="cnv-B",
        feature_label="snp",
        suffix=f"_{assay_type}",
        run_id=run_id,
        sample_id=sample_id,
        pdf=pdf,
    )
    plot_allele_freqs(
        bb_df,
        dataset_ids,
        tot_mtx_bb,
        b_mtx_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_rep_ids=cell_rep_ids,
        allele="cnv-B",
        feature_label="bb",
        suffix=f"_{assay_type}",
        run_id=run_id,
        sample_id=sample_id,
        pdf=pdf,
    )
logging.info(f"saved 2-page BAF PDF to {pdf_path}")

if is_rna_assay:
    adata: sc.AnnData = sc.read_h5ad(h5ad_file)
    barcodes = np.asarray(read_barcodes(all_barcodes), dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    assert len(missing) == 0, (
        f"h5ad, {len(missing)} barcode(s) missing, e.g. {missing[:5]}"
    )
    adata = adata[barcodes, :].copy()

    adata = assign_features_to_ranges(adata, bb_df, assay_type, range_id="bb_id")
    counts = adata.var["bb_id"].value_counts()
    bb_df["#feature"] = bb_df["bb_id"].map(counts).fillna(0).astype(int)
    x_count = sum_features_to_bbs(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    logging.info(
        f"gene-level matrix: shape={adata.X.shape}, sparsity={_sparsity(adata.X):.4f}"
    )
else:
    # scATAC: per-cell Xcount from raw 10x fragments (no tile h5ad)
    bb_grid = bb_df[["#CHR", "START", "END", "bb_id"]].copy()
    x_count = sum_atac_fragments_to_bins(
        frag_files,
        dataset_ids,
        read_full_barcodes(barcodes_full_path),
        bb_grid,
        num_bbs,
    )
logging.info(
    f"bb-level X matrix: shape={x_count.shape}, sparsity={_sparsity(x_count):.4f}"
)

save_npz(out_x_count, x_count.astype(COUNT_DTYPE))
save_npz(out_tot_mtx_bb, tot_mtx_bb.astype(COUNT_DTYPE))
save_npz(out_a_mtx_bb, a_mtx_bb.astype(COUNT_DTYPE))
save_npz(out_b_mtx_bb, b_mtx_bb.astype(COUNT_DTYPE))
bb_df.to_csv(cnv_segments, header=True, sep="\t", index=False)
shutil.copy2(all_barcodes, barcodes_out)
shutil.copy2(barcodes_full_path, barcodes_full_out)
shutil.copy2(sample_file, out_sample_file)
logging.info("finished.")
