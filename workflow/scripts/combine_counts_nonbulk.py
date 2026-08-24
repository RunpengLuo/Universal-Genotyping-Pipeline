"""Single-cell: joint adaptive binning over all non-bulk assays of the sample.

Last update: 2026-08-13

Inputs:
- allele_dir/snps.tsv.gz: the union SNP set, matrix rows
- allele_dir/snp.{T,A,B}allele.npz: one matrix over every assay's cells
- allele_dir/barcodes.tsv.gz: the matrix column axis
- allele_dir/sample_ids.tsv: dataset x assay roster, one pseudobulk each
- aux_dir/windows.bed.gz: the fixed bins, shared with bulk
- atac_fragments.tsv.gz: scATAC Xcount source, counted through the windows
- bb_dir/{assay}.h5ad: RNA Xcount source, genes assigned to bb hulls
- phase_dir/genetic_map.tsv.gz: optional, for cM-based switch probabilities
- genome_size: chrom sizes TSV
Outputs:
- bb_dir/unit/{assay}/snp.tsv.gz: the SNPs that landed in a window, matrix rows
- bb_dir/unit/{assay}/snp.{T,A,B}allele.npz: this assay's slice of their allele counts
- bb_dir/unit/{assay}/barcodes.tsv.gz: this assay's cells, matrix column order
- bb_dir/unit/{assay}/sample_ids.tsv: this assay's datasets
- bb_dir/unit/scATAC/window.{tsv.gz,Xcount.npz}: fragments counted per window per cell
- bb_dir/unit/{rna_assay}/gene.{tsv.gz,Xcount.npz}: UMIs per gene per cell, un-binned
- bb_dir/MSR{msr}/{assay}/bb.tsv.gz: shared bb definitions, duplicated per assay
- bb_dir/MSR{msr}/{assay}/bb.{T,A,B}allele.npz: this assay's per-bb allele counts
- bb_dir/MSR{msr}/{assay}/bb.Xcount.npz: this assay's per-bb native counts
- bb_dir/MSR{msr}/{assay}/barcodes.tsv.gz: this assay's cells, matrix column order
- bb_dir/MSR{msr}/{assay}/sample_ids.tsv: this assay's datasets
- bb_dir/multi_snp/{assay}/bb.tsv.gz: multi-SNP diagnostic groups, MSR-independent
- bb_dir/multi_snp/{assay}/bb.{T,A,B}allele.npz: per-group allele counts
- qc_dir/combine_counts.{assay}.MSR{msr}.pdf: allele-frequency QC per assay
"""

import logging


snakemake_handle = snakemake

from utils import log_hist, log_ratios, maybe_path, set_omp_threads, setup_logging

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from scipy.sparse import save_npz

from const import ASSAY_TYPE2MODALITY
from io_utils import (
    read_allele_mat,
    read_barcodes_by_dataset,
    read_window_bed,
    write_bb_file,
)
from combine_counts_utils import observation_cluster_ids, tumor_observation_indices
from phasing_utils import (
    estimate_switchprobs_PS,
    estimate_switchprobs_cM,
    interp_cM_between_bbs,
)
from segmentation_utils import (
    build_adaptive_bins,
    sum_features_to_bbs,
    sum_observations_to_pseudobulk,
)
from feature_utils import (
    explode_feature_ids,
    merge_feature_ids,
    read_gene_counts,
    sum_atac_fragments_to_bins,
    sum_umis_to_bins,
)
from range_utils import assign_pos_to_range, merge_ranges_to_clusters
from plot_alleles import plot_allele_freqs
from matplotlib.backends.backend_pdf import PdfPages

##################################################

# inputs
snp_info = snakemake_handle.input["snp_info"]
tot_mtx_snp_file = snakemake_handle.input["tot_mtx_snp"]
a_mtx_snp_file = snakemake_handle.input["a_mtx_snp"]
b_mtx_snp_file = snakemake_handle.input["b_mtx_snp"]
sample_file = snakemake_handle.input["sample_file"]
barcode_file = snakemake_handle.input["all_barcodes"]
frag_files = list(snakemake_handle.input["frag_files"])
h5ad_files = list(snakemake_handle.input["h5ad_files"])
gmap_file = maybe_path(snakemake_handle.input["gmap_file"])
window_bed = snakemake_handle.input["window_bed"]
genome_size = snakemake_handle.input["genome_size"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
run_id = snakemake_handle.params["run_id"]
assay_types = list(snakemake_handle.params["assay_types"])
chroms = list(snakemake_handle.params["chroms"])
nu = float(snakemake_handle.params["nu"])
min_switchprob = float(snakemake_handle.params["min_switchprob"])
switchprob_ps = float(snakemake_handle.params["switchprob_ps"])
nsnp_multi = int(snakemake_handle.params["nsnp_multi"])
msr_list = [int(m) for m in snakemake_handle.params["min_snp_reads"]]
min_snp_per_bin = int(snakemake_handle.params["min_snp_per_bin"])
gene_aware_binning = bool(snakemake_handle.params["gene_aware_binning"])

# outputs
out_bb_file = list(snakemake_handle.output["bb_file"])
out_sample_file = list(snakemake_handle.output["sample_file"])
out_tot_mtx_bb = list(snakemake_handle.output["tot_mtx_bb"])
out_a_mtx_bb = list(snakemake_handle.output["a_mtx_bb"])
out_b_mtx_bb = list(snakemake_handle.output["b_mtx_bb"])
out_unit_snp_file = list(snakemake_handle.output["unit_snp_file"])
out_unit_tot_mtx = list(snakemake_handle.output["unit_tot_mtx"])
out_unit_a_mtx = list(snakemake_handle.output["unit_a_mtx"])
out_unit_b_mtx = list(snakemake_handle.output["unit_b_mtx"])
out_unit_barcodes = list(snakemake_handle.output["unit_barcodes"])
out_unit_sample_file = list(snakemake_handle.output["unit_sample_file"])
out_unit_window_file = list(snakemake_handle.output["unit_window_file"])
out_unit_window_x = list(snakemake_handle.output["unit_window_x"])
out_unit_gene_file = list(snakemake_handle.output["unit_gene_file"])
out_unit_gene_x = list(snakemake_handle.output["unit_gene_x"])
out_multi_snp_file = list(snakemake_handle.output["multi_snp_file"])
out_tot_mtx_multi = list(snakemake_handle.output["tot_mtx_multi"])
out_a_mtx_multi = list(snakemake_handle.output["a_mtx_multi"])
out_b_mtx_multi = list(snakemake_handle.output["b_mtx_multi"])
out_all_barcodes = list(snakemake_handle.output["all_barcodes"])
out_x_count = list(snakemake_handle.output["x_count"])
out_qc_pdf = list(snakemake_handle.output["qc_pdf"])

n_assays = len(assay_types)

##################################################
# load the union: snps.tsv.gz are the rows, barcodes.tsv.gz the columns, and
# sample_ids.tsv the (dataset_id, assay_type) roster -- one row per pseudobulk column
joint_sids = pd.read_table(sample_file)
snps = pd.read_table(snp_info, sep="\t")
tot_mtx_snp = read_allele_mat(tot_mtx_snp_file)
a_mtx_snp = read_allele_mat(a_mtx_snp_file)
b_mtx_snp = read_allele_mat(b_mtx_snp_file)
n_snps = len(snps)

cells = read_barcodes_by_dataset(barcode_file)
assert len(cells) == tot_mtx_snp.shape[1], (
    f"barcodes.tsv.gz has {len(cells)} rows for {tot_mtx_snp.shape[1]} matrix columns"
)
cell_dataset_idx = observation_cluster_ids(cells, joint_sids)
assay_cols = {at: (cells["assay_type"] == at).to_numpy() for at in assay_types}
roster_rows = {at: (joint_sids["assay_type"] == at).to_numpy() for at in assay_types}
empty = [at for at, m in assay_cols.items() if not m.any()]
assert not empty, f"barcodes.tsv.gz, no cell for assay(s) {empty}"
tumor_dataset_indices = tumor_observation_indices(joint_sids)

logging.info(
    f"combine_counts_nonbulk\n"
    f"sample_id={sample_id}\n"
    f"assay_types={assay_types}\n"
    f"#SNPs={n_snps}\n"
    f"#cells={len(cells)}\n"
    f"#cells(per assay)={[int(assay_cols[at].sum()) for at in assay_types]}\n"
    f"#datasets={len(joint_sids)}\n"
    f"#tumor_datasets={len(tumor_dataset_indices)}"
)
logging.info(f"gene_aware_binning={gene_aware_binning}")

##################################################
# one pseudobulk column per (dataset_id x assay_type), in the roster's row order
tot_pb = sum_observations_to_pseudobulk(tot_mtx_snp, cell_dataset_idx, len(joint_sids))
logging.info(
    f"binning on {len(tumor_dataset_indices)}/{tot_pb.shape[1]} tumor "
    "(dataset x assay) pseudobulk observations"
)

##################################################
# shared precompute (genetic map, RNA h5ads, fixed bins)
genetic_map = pd.read_table(gmap_file, sep="\t") if gmap_file is not None else None
rna_assay_types = [at for at in assay_types if ASSAY_TYPE2MODALITY[at] == "RNA"]
assert len(h5ad_files) == len(rna_assay_types), (
    f"h5ad_files, {len(h5ad_files)} files for {len(rna_assay_types)} RNA assays"
)
h5ad_by_assay = dict(zip(rna_assay_types, h5ad_files))
# fixed bins: the window BED, the same grid bulk bins on. Tiled per segment row, so a
# window never spans two segments and none lies in a blacklist hole. Windows exist where
# no SNP does, so a SNP-free segment still yields bbs carrying Xcount.
bin_df = read_window_bed(window_bed, chroms=chroms)
logging.info(f"fixed bins: {len(bin_df)} windows from {window_bed}")

##################################################
# assign SNPs to windows, then drop the misses from the matrices too
tot_tumor = np.ascontiguousarray(tot_pb[:, tumor_dataset_indices])
snps_binned, off_idx = assign_pos_to_range(snps, bin_df, ref_id="bin_id", dropna=True)
snp_spans = (snps["END"] - snps["START"]).to_numpy()
log_ratios(
    "SNPs outside every window", len(off_idx), len(snps), snp_spans[off_idx], "snp"
)
if len(off_idx):
    log_hist(tot_tumor[off_idx].sum(axis=1), "depth of SNPs outside every bin")
log_hist(
    snps_binned.groupby("bin_id").size().reindex(range(len(bin_df)), fill_value=0),
    "SNPs per window",
)
keep_snps = np.ones(len(snps), dtype=bool)
keep_snps[off_idx] = False
tot_mtx_snp, a_mtx_snp, b_mtx_snp = (
    tot_mtx_snp[keep_snps],
    a_mtx_snp[keep_snps],
    b_mtx_snp[keep_snps],
)
tot_tumor = np.ascontiguousarray(tot_tumor[keep_snps])

##################################################
# unit level: the shared SNP grid plus each assay's native count unit, before any merge.
# scATAC counts fragments through the windows; RNA keeps whole genes, the finest grid a
# gene can be assigned to.
atac_assays = [at for at in assay_types if at == "scATAC"]
rna_assays = [at for at in assay_types if ASSAY_TYPE2MODALITY[at] == "RNA"]
unit_snps = snps_binned.drop(columns=["bin_id"])
unit_windows = bin_df.drop(columns=["bin_id"])
window_ranges = bin_df[["#CHR", "START", "END", "bin_id"]].rename(
    columns={"bin_id": "bb_id"}
)
for k, assay in enumerate(assay_types):
    cols = assay_cols[assay]
    at_sids = joint_sids[roster_rows[assay]]
    at_cells = cells[cols]
    unit_snps.to_csv(out_unit_snp_file[k], sep="\t", index=False)
    save_npz(out_unit_tot_mtx[k], tot_mtx_snp[:, cols])
    save_npz(out_unit_a_mtx[k], a_mtx_snp[:, cols])
    save_npz(out_unit_b_mtx[k], b_mtx_snp[:, cols])
    at_cells["BARCODE"].to_csv(
        out_unit_barcodes[k], sep="\t", header=False, index=False
    )
    at_sids.to_csv(out_unit_sample_file[k], sep="\t", index=False)
    logging.info(
        f"{assay}: unit level, {len(unit_snps)} SNPs x {len(at_cells)} cells to "
        f"{out_unit_snp_file[k]}"
    )
    if assay == "scATAC":
        j = atac_assays.index(assay)
        assert len(frag_files) == len(at_sids), (
            f"frag_files, {len(frag_files)} files for {len(at_sids)} scATAC datasets"
        )
        unit_windows.to_csv(out_unit_window_file[j], sep="\t", index=False)
        window_x = sum_atac_fragments_to_bins(
            frag_files,
            at_sids["dataset_id"].tolist(),
            at_cells,
            window_ranges,
            len(bin_df),
        )
        save_npz(out_unit_window_x[j], window_x)
        logging.info(
            f"{assay}: unit Xcount (fragments): shape={window_x.shape}, nnz={window_x.nnz}"
        )
    elif assay in h5ad_by_assay:
        j = rna_assays.index(assay)
        genes, gene_x = read_gene_counts(
            h5ad_by_assay[assay], at_cells["BARCODE"].tolist()
        )
        genes.to_csv(out_unit_gene_file[j], sep="\t", index=False)
        save_npz(out_unit_gene_x[j], gene_x)
        logging.info(
            f"{assay}: unit Xcount (genes): shape={gene_x.shape}, nnz={gene_x.nnz}"
        )

##################################################
# adaptive segmentation bounderies
cluster_cols = ["region_id", "seg_id"]

if "PS" in snps.columns:
    assert snps["PS"].notna().all(), "SNP file, `PS` column has NaNs"
    cluster_cols.append("PS")
    modal = snps_binned.groupby("bin_id")["PS"].agg(lambda x: x.mode().iloc[0])
    bin_df["PS"] = bin_df["bin_id"].map(modal).ffill().bfill().fillna(1)

if gene_aware_binning:
    gene_spans = (
        explode_feature_ids(snps_binned, cols=["bin_id"])
        .groupby("feature_id")["bin_id"]
        .agg(["min", "max"])
    )
    bin_df["gene_cluster"] = merge_ranges_to_clusters(
        len(bin_df), zip(gene_spans["min"].to_numpy(), gene_spans["max"].to_numpy() + 1)
    )
    logging.info(
        f"gene-aware binning: {len(gene_spans)} genes over {len(bin_df)} fixed bins -> "
        f"{bin_df['gene_cluster'].nunique()} clusters"
    )

##################################################
# multi-SNP groups: nsnp_multi SNPs each, independent of the binning sweep. One SNP per
# bin and a zero read threshold, so the grouping follows the SNP count alone and is the
# same for every assay; only the column slice of the matrices differs.
multi_cols = ["#CHR", "START", "END", "region_id"] + (
    ["seg_id"] if "seg_id" in snps_binned.columns else []
)
multi_bins = snps_binned[multi_cols].reset_index(drop=True)
multi_bins["bin_id"] = np.arange(len(multi_bins))
multi_snps_in, off_multi = assign_pos_to_range(
    snps_binned, multi_bins, ref_id="bin_id", dropna=True
)
keep_multi = np.ones(len(multi_bins), dtype=bool)
keep_multi[off_multi] = False
multi_bin_spans = (multi_bins["END"] - multi_bins["START"]).to_numpy()
log_ratios(
    "SNPs with no per-SNP range, dropped from the multi-SNP grouping",
    len(off_multi),
    len(multi_bins),
    multi_bin_spans[off_multi],
    "snp",
)
if len(off_multi):
    log_hist(tot_tumor[off_multi].sum(axis=1), "depth of SNPs with no per-SNP range")

multi_bbs, snps_multi = build_adaptive_bins(
    multi_bins,
    multi_snps_in,
    np.ascontiguousarray(tot_tumor[keep_multi]),
    0,
    nsnp_multi,
    cluster_cols=[c for c in ("region_id", "seg_id") if c in multi_bins.columns],
    max_blocksize=0,
    gene_aware=False,
)
num_multi = len(multi_bbs)
multi_ids = snps_multi["bb_id"].to_numpy()

if genetic_map is not None:
    multi_bbs["switchprobs"] = estimate_switchprobs_cM(
        interp_cM_between_bbs(multi_bbs, snps_multi, genetic_map, bb_id_col="bb_id"),
        nu=nu,
        min_switchprob=min_switchprob,
    )
else:
    multi_bbs["switchprobs"] = estimate_switchprobs_PS(multi_bbs, switchprob_ps)
multi_bbs["feature_id"] = (
    multi_bbs["bb_id"]
    .map(snps_multi.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
    .fillna("intergenic")
)

multi_cache = []
for k in range(n_assays):
    cols = assay_cols[assay_types[k]]
    tot_multi = sum_features_to_bbs(
        tot_mtx_snp[keep_multi][:, cols], multi_ids, num_multi
    )
    a_multi = sum_features_to_bbs(a_mtx_snp[keep_multi][:, cols], multi_ids, num_multi)
    b_multi = sum_features_to_bbs(b_mtx_snp[keep_multi][:, cols], multi_ids, num_multi)
    write_bb_file(multi_bbs, out_multi_snp_file[k])
    save_npz(out_tot_mtx_multi[k], tot_multi)
    save_npz(out_a_mtx_multi[k], a_multi)
    save_npz(out_b_mtx_multi[k], b_multi)
    logging.info(
        f"{assay_types[k]}: wrote {num_multi} multi-SNP groups over "
        f"{int(keep_multi.sum())} SNPs to {out_multi_snp_file[k]}"
    )
    multi_cache.append({"tot": tot_multi, "b": b_multi})

##################################################
# per-MSR joint segmentation + per-assay outputs
n_msr = len(msr_list)
for j, min_snp_reads in enumerate(msr_list):
    logging.info(f"===== MSR={min_snp_reads} =====")
    min_snp_reads_vec = np.full(
        len(tumor_dataset_indices), min_snp_reads, dtype=np.float64
    )
    bbs, snps_bb = build_adaptive_bins(
        bin_df,
        snps_binned,
        tot_tumor,
        min_snp_reads_vec,
        min_snp_per_bin,
        cluster_cols=cluster_cols,
        max_blocksize=0,
        gene_aware=gene_aware_binning,
    )
    num_bbs = len(bbs)
    bb_spans = bbs["BLOCKSIZE"].to_numpy()
    empty_bb = (bbs["#SNPS"] == 0).to_numpy()
    logging.info(f"{num_bbs} bbs from {len(bin_df)} bins")
    log_ratios(
        "SNP-free bbs (Xcount only, all-zero allele rows)",
        int(empty_bb.sum()),
        num_bbs,
        bb_spans[empty_bb],
        "bb",
    )

    if genetic_map is not None:
        dist_cms = interp_cM_between_bbs(bbs, snps_bb, genetic_map, bb_id_col="bb_id")
        bbs["switchprobs"] = estimate_switchprobs_cM(
            dist_cms, nu=nu, min_switchprob=min_switchprob
        )
    else:
        bbs["switchprobs"] = estimate_switchprobs_PS(bbs, switchprob_ps)

    bbs["feature_id"] = (
        bbs["bb_id"]
        .map(snps_bb.groupby("bb_id")["feature_id"].agg(merge_feature_ids))
        .fillna("intergenic")
    )

    # the SNP grid is shared, so every assay's matrix rows are these bb_ids
    bb_ids = snps_bb["bb_id"].to_numpy()
    # ATAC: the windows, stamped with their owning bb. A fragment in a blacklist hole or
    # between segments hits no window and is dropped, where the bb hull would swallow it.
    window_bb_ranges = bin_df[["#CHR", "START", "END", "bb_id"]].copy()
    # RNA: the bb hulls. A gene is never split, and against 1 kb windows every window
    # inside a gene ties on overlap, so gene -> window would pick one arbitrarily.
    bb_ranges = bbs[["#CHR", "START", "END", "bb_id"]].copy()

    for k in range(n_assays):
        assay = assay_types[k]
        # outputs are expanded assay-major over msr_list, so this is (k, j)
        idx = k * n_msr + j
        at_sids = joint_sids[roster_rows[assay]]
        at_cells = cells[assay_cols[assay]]
        at_cell_dataset_idx = observation_cluster_ids(at_cells, at_sids)

        write_bb_file(bbs, out_bb_file[idx])
        at_sids.to_csv(out_sample_file[idx], sep="\t", index=False)

        cols = assay_cols[assay]
        tot_bb = sum_features_to_bbs(tot_mtx_snp[:, cols], bb_ids, num_bbs)
        a_bb = sum_features_to_bbs(a_mtx_snp[:, cols], bb_ids, num_bbs)
        b_bb = sum_features_to_bbs(b_mtx_snp[:, cols], bb_ids, num_bbs)
        save_npz(out_tot_mtx_bb[idx], tot_bb)
        save_npz(out_a_mtx_bb[idx], a_bb)
        save_npz(out_b_mtx_bb[idx], b_bb)

        # per-cell Xcount per bb bin: scATAC from raw fragments, RNA from the h5ad
        if assay == "scATAC":
            assert len(frag_files) == len(at_sids), (
                f"frag_files, {len(frag_files)} files for {len(at_sids)} scATAC datasets"
            )
            x_count = sum_atac_fragments_to_bins(
                frag_files,
                at_sids["dataset_id"].tolist(),
                at_cells,
                window_bb_ranges,
                num_bbs,
            )
            save_npz(out_x_count[idx], x_count)
            logging.info(
                f"{assay} MSR={min_snp_reads} Xcount (fragments): shape={x_count.shape}, nnz={x_count.nnz}"
            )
        elif assay in h5ad_by_assay:
            x_count = sum_umis_to_bins(
                h5ad_by_assay[assay],
                at_cells["BARCODE"].tolist(),
                bb_ranges,
                num_bbs,
                assay,
            )
            save_npz(out_x_count[idx], x_count)
            logging.info(
                f"{assay} MSR={min_snp_reads} Xcount (h5ad): shape={x_count.shape}, nnz={x_count.nnz}"
            )

        with PdfPages(out_qc_pdf[idx]) as pdf:
            plot_allele_freqs(
                bbs,
                at_sids["dataset_id"].tolist(),
                at_sids["assay_type"].tolist(),
                at_sids["sample_type"].tolist(),
                tot_bb,
                b_bb,
                genome_size,
                qc_dir,
                apply_pseudobulk=True,
                cell_dataset_ids=at_cell_dataset_idx,
                allele="B",
                feature_label="bb",
                run_id=f"{assay}.MSR{min_snp_reads}.{run_id}",
                name_prefix="combine_counts",
                sample_id=sample_id,
                pdf=pdf,
            )
            multi = multi_cache[k]
            plot_allele_freqs(
                multi_bbs,
                at_sids["dataset_id"].tolist(),
                at_sids["assay_type"].tolist(),
                at_sids["sample_type"].tolist(),
                multi["tot"],
                multi["b"],
                genome_size,
                qc_dir,
                apply_pseudobulk=True,
                cell_dataset_ids=at_cell_dataset_idx,
                allele="B",
                feature_label="multi-snp",
                run_id=f"{assay}.MSR{min_snp_reads}.{run_id}",
                name_prefix="combine_counts",
                sample_id=sample_id,
                pdf=pdf,
            )

        at_cells["BARCODE"].to_csv(
            out_all_barcodes[idx], sep="\t", header=False, index=False
        )

logging.info("finished combine_counts_nonbulk.")
