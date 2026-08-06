"""Joint per-sample adaptive binning for non-bulk assays.

All non-bulk assays present in the sample (e.g. scRNA + scATAC for multiome) are segmented
on ONE shared bin grid. Each (replicate x assay) is pseudobulked into one column;
``adaptive_segmentation`` then requires ``min_snp_reads`` in EVERY column, so bins jointly
satisfy every (rep, assay) -- exactly the bulk multi-sample pattern (see combine_counts.py),
but with single cells pseudobulked per replicate first.

All outputs live under ``bb_dir/MSR{msr}/{assay}/``: the shared grid ``bb.tsv.gz`` and combined
``sample_ids.tsv`` (duplicated per assay) plus the per-assay (n_bins x n_cells) matrices
(``bb.{T,A,B}allele.npz``, ``multi_snp.*``, ``barcodes*``). Input for HATCHet3 and CalicoST.
"""

import logging
import shutil


snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from scipy.sparse import save_npz, load_npz

from io_utils import read_barcodes, read_full_barcodes
from combine_counts_utils import (
    build_union_snp_grid,
    cell_rep_idx_from_mapping,
)
from phasing_utils import (
    estimate_switchprobs_PS,
    estimate_switchprobs_cM,
    interp_cM_blocks,
    setup_phaseset_groups,
)
from matrix_utils import matrix_segmentation, pseudobulk_by_groups
from atac_utils import atac_fragments_to_bb
from rna_utils import rna_h5ad_to_bb
from aggregation_utils import (
    adaptive_segmentation,
    gene_block_labels,
    merge_feature_ids,
    snps_to_windows,
)
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
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]

# parameters
frag_reps = list(snakemake_handle.params["frag_reps"])  # parallel to frag_files
h5ad_assays = list(
    snakemake_handle.params["h5ad_assays"]
)  # parallel to h5ad_files (RNA-family)
qc_dir = snakemake_handle.params["qc_dir"]
run_id = snakemake_handle.params["run_id"]
nonbulk_assays = list(snakemake_handle.params["nonbulk_assays"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])
nsnp_multi = int(snakemake_handle.params["nsnp_multi"])
msr_list = [int(m) for m in snakemake_handle.params["min_snp_reads"]]
min_snp_per_bin = int(snakemake_handle.params["min_snp_per_bin"])
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

sample_id = (
    sample_ids_list[0]["SAMPLE_NAME"].iloc[0]
    if "SAMPLE_NAME" in sample_ids_list[0]
    else ""
)
logging.info(f"joint non-bulk binning: sample_id={sample_id}, assays={nonbulk_assays}")

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
        snps["feature_id"].notna() & (snps["feature_id"] != "intergenic"),
        ["feature_id"],
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
        .merge(snps[["#CHR", "POS0", "snp_row"]], on=["#CHR", "POS0"], how="left")[
            "snp_row"
        ]
        .to_numpy()
        .astype(np.int64)
    )
    tot_pb[shared_row, offset : offset + tot_pb_k.shape[1]] = tot_pb_k
    col_assay += [nonbulk_assays[k]] * n_reps_k
    col_repid += rep_ids_list[k]
    col_offsets.append(offset)
    offset += n_reps_k
logging.info(f"binning on {total_cols} (replicate x assay) pseudobulk columns")

##################################################
# 3. shared precompute (genetic map, sample sheet, windows, multi-SNP groups)
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None
h5ad_of = dict(zip(h5ad_assays, h5ad_files))
# combined sample sheet: one row per (replicate x assay) column
joint_sids = pd.concat(
    [sample_ids_list[k].assign(assay_type=nonbulk_assays[k]) for k in range(n_assays)],
    ignore_index=True,
)
# per-SNP zero-width windows for the joint segmentation
win_cols = ["#CHR", "START", "END", "region_id", "PS"]
if gene_aware_binning:
    win_cols.append("gene_block")
snp_windows = snps[win_cols].copy()
snp_windows["win_idx"] = np.arange(len(snp_windows))
tot_pb_cont = np.ascontiguousarray(tot_pb)
# assigned once here; every MSR below reuses it
snps_win = snps_to_windows(snps, snp_windows, tot_pb_cont)

# per-assay multi-SNP pre-grouping (diagnostic; every nsnp_multi SNPs)
multi_cache = []
for k in range(n_assays):
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
    if genetic_map is not None:
        dist_cms_multi = interp_cM_blocks(
            multi_snps, snps_multi, genetic_map, block_id_col="bb_id"
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
    logging.info(f"===== joint non-bulk binning MSR={min_snp_reads} =====")
    bbs, snps_bb = adaptive_segmentation(
        snp_windows,
        snps_win.copy(),
        tot_pb_cont,
        min_snp_reads,
        min_snp_per_bin,
        grp_cols=grp_cols,
        tumor_sidx=0,
        max_blocksize=0,
        gene_aware=gene_aware_binning,
    )
    num_bbs = len(bbs)

    if genetic_map is not None:
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

    # union SNP -> shared bb_id map + shared bb grid (fragment/gene -> bin assignment)
    bb_of_snp = snps_bb[["#CHR", "POS0", "bb_id"]]
    bb_grid = bbs[["#CHR", "START", "END", "bb_id"]].copy()

    for k in range(n_assays):
        assay = nonbulk_assays[k]
        idx = k * n_msr + j

        bb_out.to_csv(out_bb_file[idx], sep="\t", header=True, index=False)
        joint_sids.to_csv(out_sample_file[idx], sep="\t", index=False)

        # map this assay's SNPs to the shared bb grid
        m = snps_list[k][["#CHR", "POS0"]].merge(
            bb_of_snp, on=["#CHR", "POS0"], how="left"
        )
        keep = m["bb_id"].notna().to_numpy()
        rows_k = np.where(keep)[0]
        bb_ids_k = m.loc[keep, "bb_id"].to_numpy().astype(np.int64)

        tot_bb = matrix_segmentation(tot_mtx_snp_list[k][rows_k], bb_ids_k, num_bbs)
        a_bb = matrix_segmentation(a_mtx_snp_list[k][rows_k], bb_ids_k, num_bbs)
        b_bb = matrix_segmentation(b_mtx_snp_list[k][rows_k], bb_ids_k, num_bbs)
        save_npz(out_tot_mtx_bb[idx], tot_bb)
        save_npz(out_a_mtx_bb[idx], a_bb)
        save_npz(out_b_mtx_bb[idx], b_bb)

        # per-cell Xcount per bb bin: scATAC from raw fragments, RNA from the h5ad
        if assay == "scATAC":
            x_count = atac_fragments_to_bb(
                frag_files,
                frag_reps,
                read_full_barcodes(barcode_full_files[k]),
                bb_grid,
                num_bbs,
            )
            save_npz(out_x_count[idx], x_count)
            logging.info(
                f"{assay} MSR={min_snp_reads} Xcount (fragments): shape={x_count.shape}, nnz={x_count.nnz}"
            )
        elif assay in h5ad_of:
            x_count = rna_h5ad_to_bb(
                h5ad_of[assay], read_barcodes(barcode_files[k]), bb_grid, num_bbs, assay
            )
            save_npz(out_x_count[idx], x_count)
            logging.info(
                f"{assay} MSR={min_snp_reads} Xcount (h5ad): shape={x_count.shape}, nnz={x_count.nnz}"
            )

        pdf = PdfPages(out_qc_pdf[idx])
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
            run_id=f"{assay}.MSR{min_snp_reads}.{run_id}",
            name_prefix="combine_counts",
            sample_id=sample_id,
            pdf=pdf,
        )

        mc = multi_cache[k]
        mc["df"].to_csv(out_multi_snp_file[idx], sep="\t", header=True, index=False)
        save_npz(out_tot_mtx_multi[idx], mc["tot"])
        save_npz(out_a_mtx_multi[idx], mc["a"])
        save_npz(out_b_mtx_multi[idx], mc["b"])
        plot_allele_freqs(
            mc["df"],
            rep_ids_list[k],
            mc["tot"],
            mc["b"],
            genome_size,
            qc_dir,
            apply_pseudobulk=True,
            cell_rep_idx=cell_rep_idx_list[k],
            allele="B",
            unit="multi-snp",
            run_id=f"{assay}.MSR{min_snp_reads}.{run_id}",
            name_prefix="combine_counts",
            sample_id=sample_id,
            pdf=pdf,
        )
        pdf.close()

        shutil.copy2(barcode_files[k], out_all_barcodes[idx])
        shutil.copy2(barcode_full_files[k], out_barcodes_full[idx])

# TODO: recommend a default MSR (elbow of lag-1 dispersion vs #bins; see
# docs/combine_counts_pseudocode.md section 4) and record the pick.
logging.info("finished joint non-bulk binning (all MSR).")
