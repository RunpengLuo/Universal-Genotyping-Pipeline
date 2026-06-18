"""Joint per-sample adaptive binning for non-bulk assays.

All non-bulk assays present in the sample (e.g. scRNA + scATAC for multiome) are segmented
on ONE shared bin grid. Each (replicate x assay) is pseudobulked into one column;
``adaptive_segmentation`` then requires ``min_snp_reads`` in EVERY column, so bins jointly
satisfy every (rep, assay) -- exactly the bulk multi-sample pattern (see combine_counts.py),
but with single cells pseudobulked per replicate first.

Outputs are flat in ``bb_dir`` with the assay encoded in the filename: one shared
``bb.tsv.gz`` plus per-assay (n_bins x n_cells) matrices ``bb.{assay}.{T,A,B}allele.npz``.
Input for HATCHet3 and CalicoST.
"""

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
from plot_utils import plot_allele_freqs

from switchprobs import *

SNP_KEY = ["#CHR", "POS0"]

##################################################
setup_logging(snakemake.log[0])

snp_info_files = list(snakemake.input["snp_info"])
tot_files = list(snakemake.input["tot_mtx_snp"])
a_files = list(snakemake.input["a_mtx_snp"])
b_files = list(snakemake.input["b_mtx_snp"])
sample_files = list(snakemake.input["sample_file"])
barcode_files = list(snakemake.input["all_barcodes"])
barcode_full_files = list(snakemake.input["barcodes_full"])

gmap_file = maybe_path(snakemake.input["gmap_file"])
region_bed = snakemake.input["region_bed"]
genome_size = snakemake.input["genome_size"]
gtf_file = maybe_path(snakemake.input["gtf_file"])

qc_dir = snakemake.params["qc_dir"]
qc_prefix = "combine_counts"
os.makedirs(qc_dir, exist_ok=True)
run_id = snakemake.params["run_id"]
nonbulk_assays = list(snakemake.params["nonbulk_assays"])
n_assays = len(nonbulk_assays)

nu = float(snakemake.params["nu"])
min_switchprob = float(snakemake.params["min_switchprob"])
switchprob_ps = float(snakemake.params["switchprob_ps"])
nsnp_multi = int(snakemake.params["nsnp_multi"])
min_snp_reads = int(snakemake.params["min_snp_reads"])
min_snp_per_block = int(snakemake.params["min_snp_per_block"])

##################################################
# load per-assay inputs
sids_list = [pd.read_table(f) for f in sample_files]
snps_list = [pd.read_table(f, sep="\t") for f in snp_info_files]
tot_list = [load_npz(f) for f in tot_files]
a_list = [load_npz(f) for f in a_files]
b_list = [load_npz(f) for f in b_files]
rep_ids_list = [s["REP_ID"].tolist() for s in sids_list]
cell_rep_idx_list = [
    cell_rep_idx_from_mapping(read_full_barcodes(bc_full), rep_ids_list[k])
    for k, bc_full in enumerate(barcode_full_files)
]

sample_name = sids_list[0]["SAMPLE_NAME"].iloc[0] if "SAMPLE_NAME" in sids_list[0] else ""
logging.info(f"joint non-bulk binning: sample={sample_name}, assays={nonbulk_assays}")

##################################################
# 1. shared SNP grid (union across assays)
has_ps = all("PS" in s.columns for s in snps_list)
has_feature = all("feature_id" in s.columns for s in snps_list)
annot_cols = (
    ["#CHR", "POS", "POS0", "START", "END", "region_id"]
    + (["PS"] if has_ps else [])
    + (["feature_id"] if has_feature else [])
)
snps = (
    pd.concat([s[annot_cols] for s in snps_list], ignore_index=True)
    .drop_duplicates(SNP_KEY)
)
snps = sort_df_chr(snps, ch="#CHR", pos="POS0").reset_index(drop=True)
snps["snp_row"] = np.arange(len(snps))
n_snps = len(snps)
logging.info(f"shared SNP set (union): {n_snps} SNPs across {n_assays} assays")

grp_cols = ["region_id"]
if not has_ps:
    logging.info("PS not in SNP columns, setting PS=1 for all SNPs")
    snps["PS"] = 1
grp_cols.append("PS")
logging.info(f"#phaseset={snps['PS'].nunique()}")

gene_aware_binning = bool(snakemake.params["gene_aware_binning"]) and has_feature
logging.info(f"gene_aware_binning={gene_aware_binning}")
if gene_aware_binning:
    # gene blocks over the union SNPs (genomically ordered) so a bin never splits a gene
    _g = snps.loc[
        snps["feature_id"].notna() & (snps["feature_id"] != "intergenic"), ["feature_id"]
    ].copy()
    _g["__i"] = _g.index.to_numpy()
    _rng = _g.groupby("feature_id")["__i"].agg(["min", "max"])
    snps["gene_block"] = gene_block_labels(
        len(snps), zip(_rng["min"].to_numpy(), _rng["max"].to_numpy())
    )

##################################################
# 2. per-(replicate x assay) pseudobulk scattered onto the shared grid
total_cols = sum(len(r) for r in rep_ids_list)
tot_pb = np.zeros((n_snps, total_cols), dtype=np.float64)
tot_pb_list = []  # per-assay (n_snps_k x n_reps_k) pseudobulk, in that assay's SNP order
col_assay, col_repid, col_offsets = [], [], []
offset = 0
for k in range(n_assays):
    n_reps_k = len(rep_ids_list[k])
    tot_pb_k = pseudobulk_by_groups(tot_list[k], cell_rep_idx_list[k], n_reps_k)
    tot_pb_list.append(tot_pb_k)
    shared_row = (
        snps_list[k][SNP_KEY]
        .merge(snps[SNP_KEY + ["snp_row"]], on=SNP_KEY, how="left")["snp_row"]
        .to_numpy()
        .astype(np.int64)
    )
    scatter_counts_to_shared_snps(tot_pb, tot_pb_k, shared_row, offset)
    col_assay += [nonbulk_assays[k]] * n_reps_k
    col_repid += rep_ids_list[k]
    col_offsets.append(offset)
    offset += n_reps_k
logging.info(f"binning on {total_cols} (replicate x assay) pseudobulk columns")

##################################################
# 3. one joint segmentation on the shared grid (per-SNP zero-width windows)
win_cols = ["#CHR", "START", "END", "region_id", "PS"]
if gene_aware_binning:
    win_cols.append("gene_block")
snp_windows = snps[win_cols].copy()
snp_windows["win_idx"] = np.arange(len(snp_windows))

bbs, snps_bb = adaptive_segmentation(
    snp_windows,
    snps.copy(),
    np.ascontiguousarray(tot_pb),
    min_snp_reads,
    min_snp_per_block,
    grp_cols=grp_cols,
    tumor_sidx=0,
    max_blocksize=0,
    gene_aware=gene_aware_binning,
)
num_bbs = len(bbs)

_split = count_split_genes(snps_bb, grp_cols)
if _split is not None:
    logging.info(
        f"gene-split sanity: {_split[0]}/{_split[1]} genes have SNPs crossing a bin "
        f"boundary (gene_aware_binning={gene_aware_binning})"
    )

logging.info("estimate bin-level switchprobs")
if gmap_file is not None:
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms = interp_cM_blocks(bbs, snps_bb, genetic_map, block_id_col="bb_id")
    bbs["switchprobs"] = estimate_switchprobs_cM(
        dist_cms, nu=nu, min_switchprob=min_switchprob
    )
else:
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bbs[["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]].to_csv(
    snakemake.output["bb_file"], sep="\t", header=True, index=False
)
# combined sample sheet: one row per (replicate x assay) column
joint_sids = pd.concat(
    [sids_list[k].assign(assay_type=nonbulk_assays[k]) for k in range(n_assays)],
    ignore_index=True,
)
joint_sids.to_csv(snakemake.output["sample_file"], sep="\t", index=False)

# union SNP -> shared bb_id map (for aggregating each assay's per-cell matrices)
bb_of_snp = snps_bb[SNP_KEY + ["bb_id"]]

##################################################
# 4. per-assay outputs on the shared grid
for k in range(n_assays):
    assay = nonbulk_assays[k]
    qc_stamp = ".".join(p for p in (assay, run_id) if p)

    # map this assay's SNPs to the shared bb grid
    m = snps_list[k][SNP_KEY].merge(bb_of_snp, on=SNP_KEY, how="left")
    keep = m["bb_id"].notna().to_numpy()
    rows_k = np.where(keep)[0]
    bb_ids_k = m.loc[keep, "bb_id"].to_numpy().astype(np.int64)

    tot_bb = matrix_segmentation(tot_list[k][rows_k], bb_ids_k, num_bbs)
    a_bb = matrix_segmentation(a_list[k][rows_k], bb_ids_k, num_bbs)
    b_bb = matrix_segmentation(b_list[k][rows_k], bb_ids_k, num_bbs)
    save_npz(snakemake.output["tot_mtx_bb"][k], tot_bb)
    save_npz(snakemake.output["a_mtx_bb"][k], a_bb)
    save_npz(snakemake.output["b_mtx_bb"][k], b_bb)
    save_npz(snakemake.output["baf_mtx_bb"][k], csr_matrix((0, 0), dtype=np.float32))

    plot_allele_freqs(
        bbs,
        rep_ids_list[k],
        tot_bb,
        b_bb,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_rep_idx=cell_rep_idx_list[k],
        allele="B",
        unit="bb",
        run_id=qc_stamp,
        name_prefix=qc_prefix,
    )

    # ---- per-assay multi-SNP pre-grouping (diagnostic; every nsnp_multi SNPs) ----
    snps_k = snps_list[k].copy()
    if "PS" not in snps_k.columns:
        snps_k["PS"] = 1
    win_k = snps_k[["#CHR", "START", "END", "region_id"]].copy()
    win_k["win_idx"] = np.arange(len(win_k))
    multi_snps, snps_multi = adaptive_segmentation(
        win_k,
        snps_k,
        np.ascontiguousarray(tot_pb_list[k]),
        0,
        nsnp_multi,
        grp_cols=["region_id"],
        tumor_sidx=0,
        max_blocksize=0,
        gene_aware=False,
    )
    mo = snps_multi["_orig_idx"].to_numpy()
    multi_ids = snps_multi["bb_id"].to_numpy()
    tot_multi = matrix_segmentation(tot_list[k][mo], multi_ids, len(multi_snps))
    a_multi = matrix_segmentation(a_list[k][mo], multi_ids, len(multi_snps))
    b_multi = matrix_segmentation(b_list[k][mo], multi_ids, len(multi_snps))
    if gmap_file is not None:
        dist_cms_multi = interp_cM_blocks(
            multi_snps, snps_multi, genetic_map, block_id_col="bb_id"
        )
        multi_snps["switchprobs"] = estimate_switchprobs_cM(
            dist_cms_multi, nu=nu, min_switchprob=min_switchprob
        )
    else:
        multi_snps["switchprobs"] = estimate_switchprobs_PS(multi_snps, switchprob_ps)
    multi_snps.rename(columns={"bb_id": "multi_id"}).to_csv(
        snakemake.output["multi_snp_file"][k], sep="\t", header=True, index=False
    )
    save_npz(snakemake.output["tot_mtx_multi"][k], tot_multi)
    save_npz(snakemake.output["a_mtx_multi"][k], a_multi)
    save_npz(snakemake.output["b_mtx_multi"][k], b_multi)
    plot_allele_freqs(
        multi_snps,
        rep_ids_list[k],
        tot_multi,
        b_multi,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        cell_rep_idx=cell_rep_idx_list[k],
        allele="B",
        unit="multi-snp",
        run_id=qc_stamp,
        name_prefix=qc_prefix,
    )

    shutil.copy2(barcode_files[k], snakemake.output["all_barcodes"][k])
    shutil.copy2(barcode_full_files[k], snakemake.output["barcodes_full"][k])

logging.info("finished joint non-bulk binning.")
