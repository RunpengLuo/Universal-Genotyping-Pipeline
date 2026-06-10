import os, logging, shutil


t = int(getattr(snakemake, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd
from scipy.sparse import save_npz, load_npz, csr_matrix

from utils import *
from io_utils import *
from aggregation_utils import *
from combine_counts_utils import *
from count_reads_utils import *
from plot_utils import plot_allele_freqs, plot_1d_sample, plot_1d_multi_sample

from switchprobs import *

##################################################
"""
Adaptive binning by MSR and MSPB for non-bulk assays.
Input for HATCHet3 and CalicoST.
"""
setup_logging(snakemake.log[0])

snp_info = snakemake.input["snp_info"]
tot_mtx_snp = snakemake.input["tot_mtx_snp"]
a_mtx_snp = snakemake.input["a_mtx_snp"]
b_mtx_snp = snakemake.input["b_mtx_snp"]

gmap_file = maybe_path(snakemake.input["gmap_file"])
all_barcodes = maybe_path(snakemake.input["all_barcodes"])
barcodes_full_path = maybe_path(snakemake.input["barcodes_full"])
qc_dir = snakemake.params["qc_dir"]
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(snakemake.params, "run_id", "")

region_bed = snakemake.input["region_bed"]
genome_size = snakemake.input["genome_size"]
gtf_file = maybe_path(snakemake.input["gtf_file"])

sample_df = pd.read_table(snakemake.input["sample_file"])
sample_name = sample_df["SAMPLE"].iloc[0]
rep_ids = sample_df["REP_ID"].tolist()
assay_type = snakemake.params["assay_type"]

cell_rep_idx = (
    cell_rep_idx_from_mapping(read_full_barcodes(barcodes_full_path), rep_ids)
    if barcodes_full_path is not None
    else None
)

##################################################
logging.info(f"adaptive binning, sample name={sample_name}, assay_type={assay_type}")
logging.info(f"rep_ids={rep_ids}")
snps = pd.read_table(snp_info, sep="\t")

tot_mtx = load_npz(tot_mtx_snp)
a_mtx = load_npz(a_mtx_snp)
b_mtx = load_npz(b_mtx_snp)
(nsnps, nsamples) = tot_mtx.shape

assert "region_id" in snps.columns, "invalid SNP file"
grp_cols = ["region_id"]
if "PS" not in snps.columns:
    logging.info("PS not in SNP columns, setting PS=1 for all SNPs")
    snps["PS"] = 1
else:
    logging.info("PS is provided")
num_phaseset = snps["PS"].nunique()
logging.info(f"#phaseset={num_phaseset}")
grp_cols.append("PS")

multi_snps = adaptive_binning(
    snps,
    0,
    int(snakemake.params["nsnp_multi"]),
    tot_mtx,
    ["region_id"],
    colname="multi_id",
)
multi_ids = snps["multi_id"].to_numpy()
a_mtx_multi = matrix_segmentation(a_mtx, multi_ids, len(multi_snps))
b_mtx_multi = matrix_segmentation(b_mtx, multi_ids, len(multi_snps))
tot_mtx_multi = matrix_segmentation(tot_mtx, multi_ids, len(multi_snps))

plot_allele_freqs(
    multi_snps,
    rep_ids,
    tot_mtx_multi,
    b_mtx_multi,
    genome_size,
    qc_dir,
    apply_pseudobulk=True,
    cell_rep_idx=cell_rep_idx,
    allele="B",
    unit="multi-snp",
    run_id=run_id,
)
multi_snps["multi_id"] = np.arange(len(multi_snps))
if gmap_file is not None:
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms_multi = interp_cM_blocks(
        multi_snps, snps, genetic_map, block_id_col="multi_id"
    )
    multi_snps["switchprobs"] = estimate_switchprobs_cM(
        dist_cms_multi,
        nu=float(snakemake.params["nu"]),
        min_switchprob=float(snakemake.params["min_switchprob"]),
    )
else:
    switchprob_ps = float(snakemake.params["switchprob_ps"])
    multi_snps["switchprobs"] = estimate_switchprobs_PS(multi_snps, switchprob_ps)
multi_snps.to_csv(snakemake.output["multi_snp_file"], sep="\t", header=True, index=False)
save_npz(snakemake.output["tot_mtx_multi"], tot_mtx_multi)
save_npz(snakemake.output["a_mtx_multi"], a_mtx_multi)
save_npz(snakemake.output["b_mtx_multi"], b_mtx_multi)

bbs = adaptive_binning(
    snps,
    int(snakemake.params["min_snp_reads"]),
    int(snakemake.params["min_snp_per_block"]),
    tot_mtx,
    grp_cols,
    colname="bb_id",
)

bb_ids = snps["bb_id"].to_numpy()
num_bbs = len(bbs)
a_mtx_bb = matrix_segmentation(a_mtx, bb_ids, num_bbs)
b_mtx_bb = matrix_segmentation(b_mtx, bb_ids, num_bbs)
tot_mtx_bb = matrix_segmentation(tot_mtx, bb_ids, num_bbs)

plot_allele_freqs(
    bbs,
    rep_ids,
    tot_mtx_bb,
    b_mtx_bb,
    genome_size,
    qc_dir,
    apply_pseudobulk=True,
    cell_rep_idx=cell_rep_idx,
    allele="B",
    unit="bb",
    run_id=run_id,
)

bbs["bb_id"] = np.arange(len(bbs))

logging.info("estimate bin-level switchprobs")
if gmap_file is not None:
    dist_cms = interp_cM_blocks(bbs, snps, genetic_map, block_id_col="bb_id")
    bbs["switchprobs"] = estimate_switchprobs_cM(
        dist_cms,
        nu=float(snakemake.params["nu"]),
        min_switchprob=float(snakemake.params["min_switchprob"]),
    )
else:
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bbs[["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]].to_csv(
    snakemake.output["bb_file"], sep="\t", header=True, index=False
)
save_npz(snakemake.output["tot_mtx_bb"], tot_mtx_bb)
save_npz(snakemake.output["a_mtx_bb"], a_mtx_bb)
save_npz(snakemake.output["b_mtx_bb"], b_mtx_bb)
save_npz(snakemake.output["baf_mtx_bb"], csr_matrix((0, 0), dtype=np.float32))
if all_barcodes is not None:
    bb_out_dir = os.path.dirname(snakemake.output["bb_file"])
    shutil.copy2(all_barcodes, os.path.join(bb_out_dir, "barcodes.tsv.gz"))
    if barcodes_full_path is not None:
        shutil.copy2(
            barcodes_full_path, os.path.join(bb_out_dir, "barcodes.full.tsv.gz")
        )
shutil.copy2(snakemake.input["sample_file"], snakemake.output["sample_file"])
logging.info("finished.")
