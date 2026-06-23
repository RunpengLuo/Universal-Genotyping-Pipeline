"""Aggregate allele/feature-level matrix into BB block-level matrix.

Input for Copy-typing.
"""

import os
import logging
import shutil


snakemake_handle = snakemake

t = int(getattr(snakemake_handle, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd
from scipy.sparse import save_npz, load_npz
import scanpy as sc
from scipy.sparse import issparse
from utils import *
from io_utils import *
from aggregation_utils import *
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_allele_freqs

COUNT_DTYPE = np.int32


def _sparsity(X):
    """Return fraction of zero entries."""
    size = X.shape[0] * X.shape[1]
    if size == 0:
        return 0.0
    nnz = X.nnz if issparse(X) else np.count_nonzero(X)
    return 1.0 - nnz / size


log_file = snakemake_handle.log[0]
setup_logging(log_file)

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
h5ad_file = snakemake_handle.input["h5ad_file"]
ranger_dirs = snakemake_handle.input["ranger_dirs"]
all_barcodes = snakemake_handle.input["all_barcodes"]
barcodes_full_path = snakemake_handle.input["barcodes_full"]
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]
sample_file = snakemake_handle.input["sample_file"]
bb_file = snakemake_handle.input["bb_file"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_name = snakemake_handle.params["sample_name"]
assay_type = snakemake_handle.params["assay_type"]
feature_type = snakemake_handle.params["feature_type"]
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
rep_ids = sample_df["REP_ID"].tolist()

is_rna_assay = ASSAY_TYPE2MODALITY[assay_type] == "RNA"
assert assay_type not in BULK_ASSAYS, "bulk sample CNV segmentation unsupported yet"

cell_rep_idx = cell_rep_idx_from_mapping(
    read_full_barcodes(barcodes_full_path), rep_ids
)

logging.info(f"cnv segmentation, sample name={sample_name}, assay_type={assay_type}")
logging.info(f"rep_ids={rep_ids}")
snps = pd.read_table(snp_info, sep="\t")

tot_mtx = load_npz(tot_mtx_snp)
a_mtx = load_npz(a_mtx_snp)
b_mtx = load_npz(b_mtx_snp)

bb_df = pd.read_table(bb_file, sep="\t")
bb_df = sort_df_chr(bb_df, pos="START")
bb_df["bb_id"] = np.arange(len(bb_df))
num_bbs = len(bb_df)
logging.info(f"#BB blocks={num_bbs}")

snps["RAW_SNP_IDX"] = np.arange(len(snps))
snps = snp_to_region(snps, bb_df, assay_type, region_id="bb_id")

if "feature_id" in snps.columns:
    bb_df["feature_id"] = (
        bb_df["bb_id"]
        .map(snps.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )

raw_snp_ids = snps["RAW_SNP_IDX"].to_numpy()
tot_mtx = tot_mtx[raw_snp_ids, :]
a_mtx = a_mtx[raw_snp_ids, :]
b_mtx = b_mtx[raw_snp_ids, :]

logging.info(
    f"SNP-level matrices: shape={tot_mtx.shape}, "
    f"tot sparsity={_sparsity(tot_mtx):.4f}, "
    f"A sparsity={_sparsity(a_mtx):.4f}, "
    f"B sparsity={_sparsity(b_mtx):.4f}"
)

bb_ids = snps["bb_id"].to_numpy()
tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, num_bbs)
a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, num_bbs)
b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, num_bbs)
assert tot_mtx_bb.shape[0] == num_bbs

logging.info(
    f"BB-level matrices: shape={tot_mtx_bb.shape}, "
    f"T sparsity={_sparsity(tot_mtx_bb):.4f}, "
    f"A sparsity={_sparsity(a_mtx_bb):.4f}, "
    f"B sparsity={_sparsity(b_mtx_bb):.4f}"
)

pdf_path = os.path.join(qc_dir, f"cnv_segmentation.af_cnv-B_{assay_type}.{assay_type}.{run_id}.pdf")
with PdfPages(pdf_path) as pdf:
    plot_allele_freqs(
        snps,
        rep_ids,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_rep_idx=cell_rep_idx,
        allele="cnv-B",
        unit="snp",
        suffix=f"_{assay_type}",
        run_id=run_id,
        pdf=pdf,
    )
    plot_allele_freqs(
        bb_df,
        rep_ids,
        tot_mtx_bb,
        b_mtx_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_rep_idx=cell_rep_idx,
        allele="cnv-B",
        unit="bb",
        suffix=f"_{assay_type}",
        run_id=run_id,
        pdf=pdf,
    )
logging.info(f"saved 2-page BAF PDF to {pdf_path}")

if is_rna_assay:
    adata: sc.AnnData = sc.read_h5ad(h5ad_file)
    barcodes = np.asarray(read_barcodes(all_barcodes), dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    assert len(missing) == 0, (
        f"Missing {len(missing)} barcodes, e.g. {missing[:5]}, bug!"
    )
    adata = adata[barcodes, :].copy()

    adata = feature_to_blocks(
        adata, bb_df, assay_type, block_idx="bb_id", drop_cols=False
    )
    counts = adata.var["bb_id"].value_counts()
    bb_df[f"#{feature_type}"] = bb_df["bb_id"].map(counts).fillna(0).astype(int)
    x_count = matrix_segmentation(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    logging.info(
        f"Feature-level matrix: shape={adata.X.shape}, sparsity={_sparsity(adata.X):.4f}"
    )
else:
    # scATAC: per-cell Xcount from raw 10x fragments (no tile h5ad)
    frag_files = [locate_atac_fragment_file(d) for d in ranger_dirs]
    bb_grid = bb_df[["#CHR", "START", "END", "bb_id"]].copy()
    x_count = atac_fragments_to_bb(
        frag_files,
        rep_ids,
        read_full_barcodes(barcodes_full_path),
        bb_grid,
        num_bbs,
    )
logging.info(
    f"BB-level X matrix: shape={x_count.shape}, sparsity={_sparsity(x_count):.4f}"
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
