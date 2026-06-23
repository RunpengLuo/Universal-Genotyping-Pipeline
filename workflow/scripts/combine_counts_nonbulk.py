"""Joint per-sample adaptive binning for non-bulk assays.

All non-bulk assays present in the sample (e.g. scRNA + scATAC for multiome) are segmented
on ONE shared bin grid. Each (replicate x assay) is pseudobulked into one column;
``adaptive_segmentation`` then requires ``min_snp_reads`` in EVERY column, so bins jointly
satisfy every (rep, assay) -- exactly the bulk multi-sample pattern (see combine_counts.py),
but with single cells pseudobulked per replicate first.

The shared grid ``bb.tsv.gz`` and combined ``sample_ids.tsv`` are flat in ``bb_dir``;
per-assay (n_bins x n_cells) matrices live under ``bb_dir/{assay}/`` (``bb.{T,A,B}allele.npz``,
``multi_snp.*``, ``barcodes*``). Input for HATCHet3 and CalicoST.
"""

import os, logging, shutil


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

from utils import *
from io_utils import *
from aggregation_utils import *
from combine_counts_utils import *
from plot_utils import plot_allele_freqs
from matplotlib.backends.backend_pdf import PdfPages

from switchprobs import *

##################################################
log_file = snakemake_handle.log[0]
setup_logging(log_file)

# inputs
snp_info_files = list(snakemake_handle.input["snp_info"])
tot_mtx_snp_files = list(snakemake_handle.input["tot_mtx_snp"])
a_mtx_snp_files = list(snakemake_handle.input["a_mtx_snp"])
b_mtx_snp_files = list(snakemake_handle.input["b_mtx_snp"])
sample_files = list(snakemake_handle.input["sample_file"])
barcode_files = list(snakemake_handle.input["all_barcodes"])
barcode_full_files = list(snakemake_handle.input["barcodes_full"])
ranger_dirs = list(snakemake_handle.input["ranger_dirs"])
h5ad_files = list(snakemake_handle.input["h5ad_files"])
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]

# parameters
ranger_assays = list(snakemake_handle.params["ranger_assays"])  # parallel to ranger_dirs
ranger_reps = list(snakemake_handle.params["ranger_reps"])
h5ad_assays = list(snakemake_handle.params["h5ad_assays"])  # parallel to h5ad_files (RNA-family)
qc_dir = snakemake_handle.params["qc_dir"]
run_id = snakemake_handle.params["run_id"]
nonbulk_assays = list(snakemake_handle.params["nonbulk_assays"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])
nsnp_multi = int(snakemake_handle.params["nsnp_multi"])
min_snp_reads = int(snakemake_handle.params["min_snp_reads"])
min_snp_per_block = int(snakemake_handle.params["min_snp_per_block"])
gene_aware_binning_param = bool(snakemake_handle.params["gene_aware_binning"])

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

n_assays = len(nonbulk_assays)

##################################################
# load per-assay inputs
sample_ids_list = [pd.read_table(f) for f in sample_files]
snps_list = [pd.read_table(f, sep="\t") for f in snp_info_files]
tot_mtx_snp_list = [load_npz(f) for f in tot_mtx_snp_files]
a_mtx_snp_list = [load_npz(f) for f in a_mtx_snp_files]
b_mtx_snp_list = [load_npz(f) for f in b_mtx_snp_files]
rep_ids_list = [s["REP_ID"].tolist() for s in sample_ids_list]
cell_rep_idx_list = [
    cell_rep_idx_from_mapping(read_full_barcodes(bc_full), rep_ids_list[k])
    for k, bc_full in enumerate(barcode_full_files)
]

sample_name = sample_ids_list[0]["SAMPLE_NAME"].iloc[0] if "SAMPLE_NAME" in sample_ids_list[0] else ""
logging.info(f"joint non-bulk binning: sample={sample_name}, assays={nonbulk_assays}")

##################################################
# 1. shared SNP grid (union across assays)
snps, has_ps, has_feature = build_union_snp_grid(snps_list)
n_snps = len(snps)
logging.info(f"shared SNP set (union): {n_snps} SNPs across {n_assays} assays")

grp_cols = setup_phaseset_groups(snps)

gene_aware_binning = gene_aware_binning_param and has_feature
logging.info(f"gene_aware_binning={gene_aware_binning}")
if gene_aware_binning:
    # gene blocks over the union SNPs (genomically ordered) so a bin never splits a
    # gene; explode the ;-joined multi-gene feature_id so each gene gets its own span
    _g = snps.loc[
        snps["feature_id"].notna() & (snps["feature_id"] != "intergenic"), ["feature_id"]
    ].copy()
    _g["__i"] = _g.index.to_numpy()
    _g["feature_id"] = _g["feature_id"].str.split(";")
    _g = _g.explode("feature_id")
    _g = _g[_g["feature_id"] != "intergenic"]
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
    tot_pb_k = pseudobulk_by_groups(tot_mtx_snp_list[k], cell_rep_idx_list[k], n_reps_k)
    tot_pb_list.append(tot_pb_k)
    shared_row = (
        snps_list[k][["#CHR", "POS0"]]
        .merge(snps[["#CHR", "POS0", "snp_row"]], on=["#CHR", "POS0"], how="left")["snp_row"]
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

count_split_genes(snps_bb, grp_cols, gene_aware_binning)

logging.info("estimate bin-level switchprobs")
if gmap_file is not None:
    genetic_map = pd.read_table(gmap_file, sep="\t")
    dist_cms = interp_cM_blocks(bbs, snps_bb, genetic_map, block_id_col="bb_id")
    bbs["switchprobs"] = estimate_switchprobs_cM(
        dist_cms, nu=nu, min_switchprob=min_switchprob
    )
else:
    bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

bb_cols = ["#CHR", "START", "END", "#SNPS", "region_id", "switchprobs"]
if "feature_id" in snps_bb.columns:
    bbs["feature_id"] = (
        bbs["bb_id"]
        .map(snps_bb.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )
    bb_cols.append("feature_id")
bb_out = bbs[bb_cols]
for bb_path in out_bb_file:
    bb_out.to_csv(bb_path, sep="\t", header=True, index=False)
# combined sample sheet: one row per (replicate x assay) column
joint_sids = pd.concat(
    [sample_ids_list[k].assign(assay_type=nonbulk_assays[k]) for k in range(n_assays)],
    ignore_index=True,
)
for sid_path in out_sample_file:
    joint_sids.to_csv(sid_path, sep="\t", index=False)

# union SNP -> shared bb_id map (for aggregating each assay's per-cell matrices)
bb_of_snp = snps_bb[["#CHR", "POS0", "bb_id"]]
# shared bb grid (for ATAC fragment / RNA gene -> bin assignment) + per-assay Xcount paths
bb_grid = bbs[["#CHR", "START", "END", "bb_id"]].copy()
xcount_out = dict(zip(nonbulk_assays, out_x_count))
h5ad_of = dict(zip(h5ad_assays, h5ad_files))

##################################################
# 4. per-assay outputs on the shared grid
for k in range(n_assays):
    assay = nonbulk_assays[k]

    # map this assay's SNPs to the shared bb grid
    m = snps_list[k][["#CHR", "POS0"]].merge(bb_of_snp, on=["#CHR", "POS0"], how="left")
    keep = m["bb_id"].notna().to_numpy()
    rows_k = np.where(keep)[0]
    bb_ids_k = m.loc[keep, "bb_id"].to_numpy().astype(np.int64)

    tot_bb = matrix_segmentation(tot_mtx_snp_list[k][rows_k], bb_ids_k, num_bbs)
    a_bb = matrix_segmentation(a_mtx_snp_list[k][rows_k], bb_ids_k, num_bbs)
    b_bb = matrix_segmentation(b_mtx_snp_list[k][rows_k], bb_ids_k, num_bbs)
    save_npz(out_tot_mtx_bb[k], tot_bb)
    save_npz(out_a_mtx_bb[k], a_bb)
    save_npz(out_b_mtx_bb[k], b_bb)

    # per-cell Xcount per bb bin: scATAC from raw fragments, RNA from the h5ad
    if assay == "scATAC":
        _idx = [i for i, a in enumerate(ranger_assays) if a == assay]
        frag_files = [locate_atac_fragment_file(ranger_dirs[i]) for i in _idx]
        x_count = atac_fragments_to_bb(
            frag_files,
            [ranger_reps[i] for i in _idx],
            read_full_barcodes(barcode_full_files[k]),
            bb_grid,
            num_bbs,
        )
        save_npz(xcount_out[assay], x_count)
        logging.info(f"{assay} Xcount (fragments): shape={x_count.shape}, nnz={x_count.nnz}")
    elif assay in h5ad_of:
        x_count = rna_h5ad_to_bb(
            h5ad_of[assay], read_barcodes(barcode_files[k]), bb_grid, num_bbs, assay
        )
        save_npz(xcount_out[assay], x_count)
        logging.info(f"{assay} Xcount (h5ad): shape={x_count.shape}, nnz={x_count.nnz}")

    pdf = PdfPages(out_qc_pdf[k])
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
        run_id=f"{assay}.{run_id}",
        name_prefix="combine_counts",
        pdf=pdf,
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
    tot_multi = matrix_segmentation(tot_mtx_snp_list[k][mo], multi_ids, len(multi_snps))
    a_multi = matrix_segmentation(a_mtx_snp_list[k][mo], multi_ids, len(multi_snps))
    b_multi = matrix_segmentation(b_mtx_snp_list[k][mo], multi_ids, len(multi_snps))
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
        out_multi_snp_file[k], sep="\t", header=True, index=False
    )
    save_npz(out_tot_mtx_multi[k], tot_multi)
    save_npz(out_a_mtx_multi[k], a_multi)
    save_npz(out_b_mtx_multi[k], b_multi)
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
        run_id=f"{assay}.{run_id}",
        name_prefix="combine_counts",
        pdf=pdf,
    )
    pdf.close()

    shutil.copy2(barcode_files[k], out_all_barcodes[k])
    shutil.copy2(barcode_full_files[k], out_barcodes_full[k])

logging.info("finished joint non-bulk binning.")
