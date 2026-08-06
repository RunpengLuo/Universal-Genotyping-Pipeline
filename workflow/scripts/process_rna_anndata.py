import os
import logging
import tempfile

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, maybe_path, sort_chroms

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import squidpy as sq

from const import RANGER_MATRIX_H5, RANGER_SPATIAL_DIR, SPATIAL_ASSAYS
from io_utils import read_BED, read_barcodes, read_genes_gtf_file
from rna_utils import feature_to_blocks

##################################################
"""
Input:
1. 10x cell/space-ranger RNA Anndata, multiple replicates
2. reference GTF file with gene_id and intervals
3. gene blacklist
4. genome size file
5. genome regions whitelist

Output:
single h5ad matrix covers all replicates with position columns
"""

# inputs
barcode_files = snakemake_handle.input["barcodes"]
matrix_h5_files = list(snakemake_handle.input["matrix_h5"])
spatial_files = list(snakemake_handle.input.get("spatial_files", []))
gtf_file = snakemake_handle.input["gtf_file"]
gene_blacklist_file = maybe_path(snakemake_handle.input["gene_blacklist_file"])
region_bed = snakemake_handle.input["region_bed"]

# parameters
assay_type = snakemake_handle.params["assay_type"]
dataset_ids = snakemake_handle.params["dataset_ids"]
# per-dataset spatial/ filenames, aligned with spatial_files
spatial_names = list(snakemake_handle.params["spatial_names"])
gene_id_colname = str(snakemake_handle.params["gene_id_colname"])
min_frac_barcodes = float(snakemake_handle.params["min_frac_barcodes"])

# outputs
out_h5ad_file = snakemake_handle.output["h5ad_file"]


def stage_ranger_dir(tmp_dir, matrix_h5, names, paths):
    """Symlink one dataset's files into a Space Ranger layout for squidpy.read.visium.

    squidpy takes a directory, while the sample file names each spatial file
    individually so remote files can be fetched. Recreate the layout it expects:
    the feature matrix at the root, the rest under spatial/. Names come from
    RANGER_* in const.py.

    Args:
        tmp_dir: Directory to populate.
        matrix_h5: Path to this dataset's feature-barcode matrix.
        names: Space Ranger filenames under spatial/, for this dataset.
        paths: Paths supplying those files, in the same order.

    Returns:
        tmp_dir, ready to pass to squidpy.read.visium.
    """
    os.symlink(
        os.path.abspath(matrix_h5),
        os.path.join(tmp_dir, RANGER_MATRIX_H5[0]),
    )
    spatial_dir = os.path.join(tmp_dir, RANGER_SPATIAL_DIR)
    os.makedirs(spatial_dir)
    for name, path in zip(names, paths):
        os.symlink(os.path.abspath(path), os.path.join(spatial_dir, name))
    return tmp_dir


logging.info(f"prepare rna anndata, assay_type={assay_type}, dataset_ids={dataset_ids}")

adatas = {}
_spatial_offset = 0
for idx, dataset_id in enumerate(dataset_ids):
    logging.info(f"process {assay_type}-{dataset_id}")
    barcodes = read_barcodes(barcode_files[idx])
    barcodes = pd.Index(barcodes).astype(str)

    matrix_h5 = matrix_h5_files[idx]
    if assay_type in SPATIAL_ASSAYS:
        names = spatial_names[idx]
        paths = spatial_files[_spatial_offset : _spatial_offset + len(names)]
        _spatial_offset += len(names)
        # squidpy doesn't support load images from 3' data yet.
        load_images = assay_type == "VISIUM"
        with tempfile.TemporaryDirectory() as tmp_dir:
            ranger_dir = stage_ranger_dir(tmp_dir, matrix_h5, names, paths)
            logging.info(f"staged {len(names) + 1} files for squidpy: {names}")
            adata: sc.AnnData = sq.read.visium(
                ranger_dir, load_images=load_images, library_id=dataset_id
            )
        adata.var_names_make_unique()
    else:
        adata: sc.AnnData = sc.read_10x_h5(matrix_h5, gex_only=True)
        adata.var_names_make_unique()

    adata.obs_names = adata.obs_names.astype(str)
    adata = adata[adata.obs_names.isin(barcodes), :].copy()
    adata.obs_names = adata.obs_names.astype(str) + f"_{dataset_id}"
    adatas[dataset_id] = adata
    logging.info(f"#barcodes={adata.n_obs}, #features={adata.n_vars}")

if len(adatas) > 1:
    adata = anndata.concat(
        adatas,
        join="outer",  # union of var (genes)
        label="REP_ID",
        merge="same",
        uns_merge="unique",
        fill_value=0,
    )
else:
    adata = adatas[dataset_ids[0]]
adata.X = adata.X.tocsr()
num_total_barcodes = adata.n_obs
logging.info(f"#concat barcodes={num_total_barcodes}, #union features={adata.n_vars}")

genes_gtf = read_genes_gtf_file(gtf_file, id_col=gene_id_colname)[
    [gene_id_colname, "#CHR", "START", "END"]
]
logging.info(f"loaded #{len(genes_gtf)} unique genes from GTF.")

var_coords = adata.var.merge(
    genes_gtf, how="left", on=gene_id_colname, validate="m:1", sort=False
)
var_coords.index = adata.var.index

na_genes = var_coords["START"].isna().to_numpy()
logging.warning(
    f"#genes not found in reference GTF file={na_genes.sum()}/{len(var_coords)}"
)

adata = adata[:, ~na_genes].copy()
adata.var = var_coords.loc[~na_genes, :].copy()

adata.var["#CHR"] = adata.var["#CHR"].astype(str)
adata.var["START"] = adata.var["START"].astype(int)
adata.var["END"] = adata.var["END"].astype(int)

adata.var["pseudobulk_counts"] = np.asarray(adata.X.sum(axis=0)).flatten()
adata = adata[:, adata.var["pseudobulk_counts"] > 0].copy()
logging.info(f"#genes after filtering by zero pseudobulk_counts: {adata.n_vars}")

if gene_blacklist_file is not None:
    gene_blacklist = (
        pd.read_table(gene_blacklist_file, header=None).iloc[:, 0].to_numpy()
    )
    ind_gene_blacklist = np.isin(adata.var.index, gene_blacklist)
    logging.info(
        f"remove #{np.sum(ind_gene_blacklist)}/{adata.n_vars} genes based on {gene_blacklist_file}"
    )
    adata = adata[:, ~ind_gene_blacklist]

sum_count_before_filtering = float(adata.X.sum())
min_expressed_barcodes = round(min_frac_barcodes * num_total_barcodes)
logging.info(
    f"min_frac_barcodes={min_frac_barcodes}, min_expressed_barcodes={min_expressed_barcodes}/{num_total_barcodes}"
)

if assay_type in SPATIAL_ASSAYS:
    nnz_per_gene = adata.X.getnnz(axis=0)
    ind_sufficient_expressed_genes = np.asarray(
        nnz_per_gene >= min_expressed_barcodes
    ).ravel()
    adata = adata[:, ind_sufficient_expressed_genes].copy()
    count_ratio = float(adata.X.sum()) / sum_count_before_filtering
    logging.info(
        f"Retaining {100.0 * np.mean(ind_sufficient_expressed_genes):.3f}% of genes with sufficient expression across spots ({100.0 * count_ratio:.2f}% of total UMIs) @ {min_frac_barcodes} fraction of barcodes."
    )

regions = read_BED(region_bed)[["#CHR", "START", "END", "region_id"]]
adata = feature_to_blocks(adata, regions, assay_type)

chs = sort_chroms(adata.var["#CHR"].unique().tolist())
adata.var["#CHR"] = pd.Categorical(adata.var["#CHR"], categories=chs, ordered=True)

assert adata.var_names.is_unique, "var_names is not unique!"
sort_index = adata.var.sort_values(by=["#CHR", "START"]).index
adata = adata[:, sort_index].copy()

adata.write_h5ad(out_h5ad_file, compression="gzip")

logging.info(f"final processed {assay_type} AnnData")
logging.info(f"final #obs={adata.n_obs}, #vars={adata.n_vars}")
