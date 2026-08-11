"""Aggregate one non-bulk assay's SNP and gene counts onto pre-computed bbs.

The `copytyping_preprocess` mode: the bbs come from `bb_file` (config), so nothing is
binned here. SNPs are assigned to their containing bb and summed; the per-cell count
matrix is the h5ad UMIs (RNA) or raw 10x fragments (scATAC). Input for Copy-typing.

Inputs
  snp_info, {tot,a,b}_mtx_snp: this assay's SNP-level matrices, from phase_and_concat
  bb_file: pre-computed bbs, the feature axis of every output
  h5ad_file / frag_files: the per-cell count source, by modality
Outputs:
  bb_file: the given bbs, re-stamped with this assay's `#SNPS` / `feature_id`, plus
    `bb_id` and (RNA) `#feature`; every other column of the input passes through
  bb.{Xcount,Tallele,Aallele,Ballele}.npz: (bb x cell) matrices
  barcodes, sample_ids: this assay's slice of the union inputs
"""

import logging


snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, sort_df_chr

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from scipy.sparse import save_npz
import scanpy as sc
from const import ASSAY_TYPE2MODALITY
from io_utils import (
    read_barcodes_by_dataset,
    read_snp_mats,
    write_bb_file,
)
from combine_counts_utils import observation_cluster_ids
from segmentation_utils import sum_features_to_bbs
from range_utils import assign_pos_to_range
from feature_utils import (
    assign_features_to_ranges,
    merge_feature_ids,
    sum_atac_fragments_to_bins,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_alleles import plot_allele_freqs

COUNT_DTYPE = np.int32

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
h5ad_file = snakemake_handle.input["h5ad_file"]
frag_files = list(snakemake_handle.input["frag_files"])
all_barcodes = snakemake_handle.input["all_barcodes"]
genome_size = snakemake_handle.input["genome_size"]
sample_file = snakemake_handle.input["sample_file"]
bb_file = snakemake_handle.input["bb_file"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
assay_type = snakemake_handle.params["assay_type"]
run_id = snakemake_handle.params["run_id"]

# outputs
out_x_count = snakemake_handle.output["x_count"]
out_tot_mtx_bb = snakemake_handle.output["tot_mtx_bb"]
out_a_mtx_bb = snakemake_handle.output["a_mtx_bb"]
out_b_mtx_bb = snakemake_handle.output["b_mtx_bb"]
out_bb_file = snakemake_handle.output["bb_file"]
out_barcodes = snakemake_handle.output["barcodes_out"]
out_sample_file = snakemake_handle.output["sample_file"]
out_qc_pdf = snakemake_handle.output["qc_pdf"]


# the allele inputs are the union over every non-bulk assay, so slice to this one
sample_df = pd.read_table(sample_file)
sample_df = sample_df[sample_df["assay_type"] == assay_type].reset_index(drop=True)
dataset_ids = sample_df["dataset_id"].tolist()
sample_types = sample_df["sample_type"].tolist()

is_rna_assay = ASSAY_TYPE2MODALITY[assay_type] == "RNA"

cells = read_barcodes_by_dataset(all_barcodes)
assay_cols = (cells["assay_type"] == assay_type).to_numpy()
assert assay_cols.any(), f"barcodes.tsv.gz, no cell for assay {assay_type}"
cells = cells[assay_cols]
cell_dataset_ids = observation_cluster_ids(cells, sample_df)

snps, tot_mtx, a_mtx, b_mtx = read_snp_mats(snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp)
tot_mtx, a_mtx, b_mtx = (m[:, assay_cols] for m in (tot_mtx, a_mtx, b_mtx))

bb_df = pd.read_table(bb_file, sep="\t")
bb_df = sort_df_chr(bb_df, pos="START")
bb_df["bb_id"] = np.arange(len(bb_df))
num_bbs = len(bb_df)

logging.info(
    f"combine_counts_fixed_bins\n"
    f"sample_id={sample_id}\n"
    f"assay_type={assay_type}\n"
    f"#SNPs={len(snps)}\n"
    f"#datasets={len(dataset_ids)}\n"
    f"#bbs={num_bbs}"
)

snps["RAW_SNP_DF_IDX"] = np.arange(len(snps))
snps, _ = assign_pos_to_range(snps, bb_df, ref_id="bb_id", dropna=True)
logging.info(f"#SNPs in a bb={len(snps)}")
bb_df["#SNPS"] = bb_df["bb_id"].map(snps["bb_id"].value_counts()).fillna(0).astype(int)

bb_df["feature_id"] = (
    bb_df["bb_id"]
    .map(snps.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
    .fillna("intergenic")
)

raw_snp_df_idx = snps["RAW_SNP_DF_IDX"].to_numpy()
tot_mtx = tot_mtx[raw_snp_df_idx, :]
a_mtx = a_mtx[raw_snp_df_idx, :]
b_mtx = b_mtx[raw_snp_df_idx, :]

logging.info(f"SNP-level matrices: shape={tot_mtx.shape}, tot nnz={tot_mtx.nnz}")

bb_ids = snps["bb_id"].to_numpy()
tot_mtx_bb = sum_features_to_bbs(tot_mtx, bb_ids, num_bbs)
a_mtx_bb = sum_features_to_bbs(a_mtx, bb_ids, num_bbs)
b_mtx_bb = sum_features_to_bbs(b_mtx, bb_ids, num_bbs)

logging.info(f"bb-level matrices: shape={tot_mtx_bb.shape}, T nnz={tot_mtx_bb.nnz}")

with PdfPages(out_qc_pdf) as pdf:
    plot_allele_freqs(
        snps,
        dataset_ids,
        [assay_type] * len(dataset_ids),
        sample_types,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_dataset_ids=cell_dataset_ids,
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
        [assay_type] * len(dataset_ids),
        sample_types,
        tot_mtx_bb,
        b_mtx_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_dataset_ids=cell_dataset_ids,
        allele="cnv-B",
        feature_label="bb",
        suffix=f"_{assay_type}",
        run_id=run_id,
        sample_id=sample_id,
        pdf=pdf,
    )
logging.info(f"saved 2-page BAF PDF to {out_qc_pdf}")

if is_rna_assay:
    adata: sc.AnnData = sc.read_h5ad(h5ad_file)
    barcodes = cells["BARCODE"].to_numpy().astype(str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    assert len(missing) == 0, (
        f"h5ad, {len(missing)} barcode(s) missing, e.g. {missing[:5]}"
    )
    adata = adata[barcodes, :].copy()

    adata = assign_features_to_ranges(adata, bb_df, assay_type, range_id="bb_id")
    counts = adata.var["bb_id"].value_counts()
    bb_df["#feature"] = bb_df["bb_id"].map(counts).fillna(0).astype(int)
    x_count = sum_features_to_bbs(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    logging.info(f"gene-level matrix: shape={adata.X.shape}")
else:
    # scATAC: per-cell Xcount from raw 10x fragments (no tile h5ad)
    bb_grid = bb_df[["#CHR", "START", "END", "bb_id"]].copy()
    x_count = sum_atac_fragments_to_bins(
        frag_files,
        dataset_ids,
        cells,
        bb_grid,
        num_bbs,
    )
logging.info(f"bb-level X matrix: shape={x_count.shape}, nnz={x_count.nnz}")

save_npz(out_x_count, x_count.astype(COUNT_DTYPE))
save_npz(out_tot_mtx_bb, tot_mtx_bb.astype(COUNT_DTYPE))
save_npz(out_a_mtx_bb, a_mtx_bb.astype(COUNT_DTYPE))
save_npz(out_b_mtx_bb, b_mtx_bb.astype(COUNT_DTYPE))
write_bb_file(bb_df, out_bb_file)
cells["BARCODE"].to_csv(out_barcodes, sep="\t", header=False, index=False)
sample_df.to_csv(out_sample_file, sep="\t", index=False)
logging.info("finished combine_counts_fixed_bins.")
