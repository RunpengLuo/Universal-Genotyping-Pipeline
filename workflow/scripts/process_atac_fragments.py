import os, logging
from snakemake.script import snakemake as sm

t = int(getattr(sm, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd
import anndata
from scanpy import AnnData
import snapatac2 as snap

from io_utils import *
from aggregation_utils import *

##################################################
"""
Input:
1. 10x cell-ranger ATAC fragments, multiple replicates
2. reference GTF file with gene_id and intervals
3. gene blacklist
4. genome size file
5. genome regions whitelist

Output:
single h5ad matrix covers all replicates with position columns
"""
setup_logging(sm.log[0])

barcode_files = sm.input["barcodes"]
ranger_dirs = sm.input["ranger_dirs"]
assay_type = sm.params["assay_type"]
rep_ids = sm.params["rep_ids"]
chrom_sizes = get_chr_sizes(sm.input["genome_size"])
tile_width = int(sm.params["tile_width"])
logging.info(f"prepare atac anndata, assay_type={assay_type}, rep_ids={rep_ids}")
logging.info(f"tile_width={tile_width}bp")

adatas = {}
for idx, rep_id in enumerate(rep_ids):
    logging.info(f"process {assay_type}-{rep_id}")
    barcodes = read_barcodes(barcode_files[idx])
    ranger_dir = ranger_dirs[idx]
    fragment_file = None
    for fpath in ["atac_fragments.tsv.gz", "fragments.tsv.gz"]:
        if os.path.exists(os.path.join(ranger_dir, fpath)):
            fragment_file = os.path.join(ranger_dir, fpath)
            break
    assert fragment_file is not None, (
        f"failed to locate atac fragment files for {rep_id}"
    )
    adata: AnnData = snap.pp.import_fragments(
        fragment_file,
        chrom_sizes,
        whitelist=barcodes,
        min_num_fragments=0,
        sorted_by_barcode=False,
        tempdir=None,
        n_jobs=int(sm.threads),
    )
    adata.obs_names = adata.obs_names.astype(str) + f"_{rep_id}"
    snap.pp.add_tile_matrix(
        adata,
        bin_size=tile_width,
        n_jobs=int(sm.threads),
    )
    adatas[rep_id] = adata
    logging.info(f"#barcodes={adata.n_obs}, #features={adata.n_vars}")

if len(adatas) > 1:
    adata = anndata.concat(
        adatas,
        join="inner",  # union of var (genes)
        label="REP_ID",
    )
else:
    adata = adatas[rep_ids[0]]

adata.X = adata.X.tocsr()
logging.info(f"#concat barcodes={adata.n_obs}, #union features={adata.n_vars}")

adata.var["pseudobulk_counts"] = np.asarray(adata.X.sum(axis=0)).flatten()
adata = adata[:, adata.var["pseudobulk_counts"] > 0].copy()

adata.var["#CHR"] = adata.var_names.str.split(":").str[0].astype(str)
chs = sort_chroms(adata.var["#CHR"].unique().tolist())
adata.var["#CHR"] = pd.Categorical(adata.var["#CHR"], categories=chs, ordered=True)
adata.var["START"] = (
    adata.var_names.str.split(":").str[1].str.split("-").str[0].astype(int)
)
adata.var["END"] = (
    adata.var_names.str.split(":").str[1].str.split("-").str[1].astype(int)
)

regions = read_region_file(sm.input["region_bed"])
adata = feature_to_blocks(adata, regions, assay_type)

assert adata.var_names.is_unique, "var_names is not unique!"
sort_index = adata.var.sort_values(by=["#CHR", "START"]).index
adata = adata[:, sort_index].copy()
# drop bulky snapatac2 fragments and metadata to reduce h5ad size
adata.obsm.clear()
adata.uns.clear()

adata.write_h5ad(sm.output["h5ad_file"], compression="gzip")

logging.info(f"final processed {assay_type} AnnData")
logging.info(f"final #obs={adata.n_obs}, #vars={adata.n_vars}")
