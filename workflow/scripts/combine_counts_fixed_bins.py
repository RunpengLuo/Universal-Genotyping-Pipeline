"""copytyping_preprocess: count one non-bulk assay onto pre-computed bbs.

Last update: 2026-08-29

Inputs:
- bb_file: pre-computed bbs, the feature axis of every output
- allele_dir/snps.tsv.gz: the union SNP set, matrix rows
- allele_dir/snp.{T,A,B}allele.npz: union allele counts, sliced to this assay
- allele_dir/barcodes.tsv.gz: union column axis, sliced to this assay
- allele_dir/sample_ids.tsv: dataset roster, sliced to this assay
- bb_dir/{assay}.h5ad: RNA UMI source for Xcount
- atac_fragments.tsv.gz: scATAC fragment source for Xcount
- aux_dir/windows.bed.gz: the window grid; scATAC fragments are counted through it
- genome_size: chrom sizes TSV
Outputs:
- bb_dir/unit/{assay}/snp.tsv.gz: the SNPs that landed in a window, matrix rows
- bb_dir/unit/{assay}/snp.{T,A,B}allele.npz: this assay's slice of their allele counts
- bb_dir/unit/{assay}/barcodes.tsv.gz: this assay's cells, matrix column order
- bb_dir/unit/{assay}/sample_ids.tsv: this assay's datasets
- bb_dir/unit/scATAC/window.{tsv.gz,Xcount.npz}: fragments counted per window per cell
- bb_dir/unit/{rna_assay}/gene.{tsv.gz,Xcount.npz}: UMIs per gene per cell, un-binned
- bb_dir/{assay}/bb.tsv.gz: the given bbs, re-stamped for this assay
- bb_dir/{assay}/bb.{Xcount,Tallele,Aallele,Ballele}.npz: per-bb count matrices
- bb_dir/{assay}/barcodes.tsv.gz: this assay's cells, matrix column order
- bb_dir/{assay}/sample_ids.tsv: this assay's datasets
- qc_dir/combine_counts_fixed_bins.{assay}.pdf: per dataset, one SNP BAF page then one
  bb page carrying pseudobulk RDR over BAF
"""

import logging


snakemake_handle = snakemake

from utils import log_ratios, set_omp_threads, setup_logging, sort_df_chr

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
    read_window_bed,
    write_bb_file,
)
from combine_counts_utils import observation_cluster_ids
from segmentation_utils import sum_features_to_bbs
from range_utils import assign_pos_to_range, assign_range_to_range
from feature_utils import (
    assign_features_to_ranges,
    merge_feature_ids,
    read_gene_counts,
    sum_atac_fragments_to_bins,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_alleles import compute_af_by_clusters
from plot_combine_counts import compute_pseudobulk_rdr, plot_pseudobulk_tracks

COUNT_DTYPE = np.int32
# the pseudobulk RDR is a share of library size, not a ratio to a reference column
RDR_YLABEL = r"$\sum_i X_{i,g}/\sum_i T_i$"

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp = snakemake_handle.input["b_mtx_snp"]
h5ad_file = snakemake_handle.input["h5ad_file"]
frag_files = list(snakemake_handle.input["frag_files"])
all_barcodes = snakemake_handle.input["all_barcodes"]
genome_size = snakemake_handle.input["genome_size"]
window_bed = snakemake_handle.input["window_bed"]
sample_file = snakemake_handle.input["sample_file"]
bb_file = snakemake_handle.input["bb_file"]

# parameters
sample_id = snakemake_handle.params["sample_id"]
assay_type = snakemake_handle.params["assay_type"]
chroms = list(snakemake_handle.params["chroms"])

# outputs
out_x_count = snakemake_handle.output["x_count"]
out_tot_mtx_bb = snakemake_handle.output["tot_mtx_bb"]
out_a_mtx_bb = snakemake_handle.output["a_mtx_bb"]
out_b_mtx_bb = snakemake_handle.output["b_mtx_bb"]
out_bb_file = snakemake_handle.output["bb_file"]
out_barcodes = snakemake_handle.output["barcodes_out"]
out_sample_file = snakemake_handle.output["sample_file"]
out_unit_snp_file = snakemake_handle.output["unit_snp_file"]
out_unit_tot_mtx = snakemake_handle.output["unit_tot_mtx"]
out_unit_a_mtx = snakemake_handle.output["unit_a_mtx"]
out_unit_b_mtx = snakemake_handle.output["unit_b_mtx"]
out_unit_barcodes = snakemake_handle.output["unit_barcodes"]
out_unit_sample_file = snakemake_handle.output["unit_sample_file"]
out_qc_pdf = snakemake_handle.output["qc_pdf"]


# the allele inputs are the union over every non-bulk assay, so slice to this one
sample_df = pd.read_table(sample_file)
sample_df = sample_df[sample_df["assay_type"] == assay_type].reset_index(drop=True)
dataset_ids = sample_df["dataset_id"].tolist()
dataset_assays = sample_df["assay_type"].tolist()
sample_types = sample_df["sample_type"].tolist()

is_rna_assay = ASSAY_TYPE2MODALITY[assay_type] == "RNA"

cells = read_barcodes_by_dataset(all_barcodes)
assay_cols = (cells["assay_type"] == assay_type).to_numpy()
assert assay_cols.any(), f"barcodes.tsv.gz, no cell for assay {assay_type}"
cells = cells[assay_cols]
cell_dataset_ids = observation_cluster_ids(cells, sample_df)

snps, tot_mtx, a_mtx, b_mtx = read_snp_mats(snp_info, tot_mtx_snp, a_mtx_snp, b_mtx_snp)
tot_mtx, a_mtx, b_mtx = (m[:, assay_cols] for m in (tot_mtx, a_mtx, b_mtx))

##################################################
# unit level: the SNP grid and this assay's native count unit, before the given bbs
bin_df = read_window_bed(window_bed, chroms=chroms)
snps_binned, off_idx = assign_pos_to_range(snps, bin_df, ref_id="bin_id", dropna=True)
snp_spans = (snps["END"] - snps["START"]).to_numpy()
log_ratios(
    "SNPs outside every window", len(off_idx), len(snps), snp_spans[off_idx], "snp"
)
keep_snps = np.ones(len(snps), dtype=bool)
keep_snps[off_idx] = False
snps_binned.drop(columns=["bin_id"]).to_csv(out_unit_snp_file, sep="\t", index=False)
save_npz(out_unit_tot_mtx, tot_mtx[keep_snps].astype(COUNT_DTYPE))
save_npz(out_unit_a_mtx, a_mtx[keep_snps].astype(COUNT_DTYPE))
save_npz(out_unit_b_mtx, b_mtx[keep_snps].astype(COUNT_DTYPE))
cells["BARCODE"].to_csv(out_unit_barcodes, sep="\t", header=False, index=False)
sample_df.to_csv(out_unit_sample_file, sep="\t", index=False)
logging.info(
    f"unit level: {len(snps_binned)} SNPs x {len(cells)} cells to {out_unit_snp_file}"
)

# the windows carry no blacklisted span, so this drops the blacklisted SNPs the given
# bb hulls would otherwise swallow
snps = snps_binned.drop(columns=["bin_id"])
tot_mtx, a_mtx, b_mtx = tot_mtx[keep_snps], a_mtx[keep_snps], b_mtx[keep_snps]

if is_rna_assay:
    genes, gene_x = read_gene_counts(h5ad_file, cells["BARCODE"].tolist())
    genes.to_csv(snakemake_handle.output["unit_gene_file"], sep="\t", index=False)
    save_npz(snakemake_handle.output["unit_gene_x"], gene_x)
    logging.info(f"unit Xcount (genes): shape={gene_x.shape}, nnz={gene_x.nnz}")
else:
    bin_df.drop(columns=["bin_id"]).to_csv(
        snakemake_handle.output["unit_window_file"], sep="\t", index=False
    )
    window_x = sum_atac_fragments_to_bins(
        frag_files,
        dataset_ids,
        cells,
        bin_df[["#CHR", "START", "END", "bin_id"]].rename(columns={"bin_id": "bb_id"}),
        len(bin_df),
    )
    save_npz(snakemake_handle.output["unit_window_x"], window_x)
    logging.info(f"unit Xcount (fragments): shape={window_x.shape}, nnz={window_x.nnz}")

##################################################
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
    # scATAC: per-cell Xcount from raw 10x fragments, counted through the windows
    # stamped with their owning bb. A fragment in a blacklist hole then hits no window
    # and is dropped, where the bb hull would swallow it.
    bb_grid, off_win = assign_range_to_range(
        bin_df[["#CHR", "START", "END"]], bb_df, "bb_id", rule="midpoint", dropna=True
    )
    win_spans = (bin_df["END"] - bin_df["START"]).to_numpy()
    log_ratios(
        "windows outside every bb", len(off_win), len(bin_df), win_spans[off_win], "bin"
    )
    x_count = sum_atac_fragments_to_bins(
        frag_files,
        dataset_ids,
        cells,
        bb_grid,
        num_bbs,
    )
logging.info(f"bb-level X matrix: shape={x_count.shape}, nnz={x_count.nnz}")

num_datasets = len(dataset_ids)
baf_snp = compute_af_by_clusters(tot_mtx, b_mtx, cell_dataset_ids, num_datasets)
baf_bb = compute_af_by_clusters(tot_mtx_bb, b_mtx_bb, cell_dataset_ids, num_datasets)
rdr_bb = compute_pseudobulk_rdr(x_count, cell_dataset_ids, num_datasets)

with PdfPages(out_qc_pdf) as pdf:
    plot_pseudobulk_tracks(
        snps,
        [("BAF", baf_snp, "BAF")],
        sample_id,
        dataset_ids,
        dataset_assays,
        sample_types,
        genome_size,
        out_qc_pdf,
        feature_label="SNP",
        pdf=pdf,
    )
    plot_pseudobulk_tracks(
        bb_df,
        [("RDR", rdr_bb, RDR_YLABEL), ("BAF", baf_bb, "BAF")],
        sample_id,
        dataset_ids,
        dataset_assays,
        sample_types,
        genome_size,
        out_qc_pdf,
        feature_label="bb",
        pdf=pdf,
    )
logging.info(f"saved QC PDF to {out_qc_pdf}")

save_npz(out_x_count, x_count.astype(COUNT_DTYPE))
save_npz(out_tot_mtx_bb, tot_mtx_bb.astype(COUNT_DTYPE))
save_npz(out_a_mtx_bb, a_mtx_bb.astype(COUNT_DTYPE))
save_npz(out_b_mtx_bb, b_mtx_bb.astype(COUNT_DTYPE))
write_bb_file(bb_df, out_bb_file)
cells["BARCODE"].to_csv(out_barcodes, sep="\t", header=False, index=False)
sample_df.to_csv(out_sample_file, sep="\t", index=False)
logging.info("finished combine_counts_fixed_bins.")
