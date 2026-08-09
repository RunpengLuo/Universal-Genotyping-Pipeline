"""Phase-and-concat allele counts for a single non-bulk assay (per assay_type).

Merges the assay's replicate pileups onto the shared parent SNP set, phases ref/alt
into A/B, filters SNPs (region, blacklist, RNA gene overlap / ATAC none), and writes
sparse (SNPs x cells) matrices plus per-cell barcode bookkeeping. Joint binning across
the sample's assays happens downstream in combine_counts_nonbulk.
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import save_npz

from const import ASSAY_TYPE2MODALITY
from io_utils import read_BED, read_VCF, write_sample_ids, write_snp_info
from range_utils import overlaps_any_range
from phase_and_concat_utils import (
    apply_masks_to_df,
    get_mask_by_blacklist,
    get_mask_by_exon,
    get_mask_by_region,
    hstack_replicate_mats,
    interp_pos_ranges,
    map_allele_mat_to_snps,
)
from matplotlib.backends.backend_pdf import PdfPages
from plot_alleles import plot_allele_freqs, plot_snp_depth
from phasing_utils import apply_phase_to_mat
from feature_utils import annotate_feature_type


##################################################
logging.info("phase and concat allele-level count matrices")

# inputs
vcf_files = snakemake_handle.input["vcfs"]
sample_tsvs = snakemake_handle.input["sample_tsvs"]
tot_mtx_files = snakemake_handle.input["tot_mtxs"]
ad_mtx_files = snakemake_handle.input["ad_mtxs"]
snp_vcf = snakemake_handle.input["snp_vcf"]
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]
gtf_file = maybe_path(snakemake_handle.input["gtf_file"])
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])
h5ad_file = snakemake_handle.input["h5ad_file"]

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
assay_type = snakemake_handle.params["assay_type"]
dataset_ids = snakemake_handle.params["dataset_ids"]
sample_types = snakemake_handle.params["sample_types"]
exon_only = snakemake_handle.params["exon_only"]
run_id = snakemake_handle.params["run_id"]

# outputs
out_snp_info = snakemake_handle.output["snp_info"]
out_tot_mtx_snp = snakemake_handle.output["tot_mtx_snp"]
out_a_mtx_snp = snakemake_handle.output["a_mtx_snp"]
out_b_mtx_snp = snakemake_handle.output["b_mtx_snp"]
out_unique_snp_ids = snakemake_handle.output["unique_snp_ids"]
out_all_barcodes = snakemake_handle.output["all_barcodes"]
out_barcodes_full = snakemake_handle.output["barcodes_full"]
out_sample_file = snakemake_handle.output["sample_file"]
out_qc_pdf = snakemake_handle.output["qc_pdf"]

is_rna_assay = ASSAY_TYPE2MODALITY[assay_type] == "RNA"

##################################################
logging.info(
    f"sample_id={sample_id}, assay_type={assay_type}, dataset_ids={dataset_ids}"
)

snps = read_VCF(snp_vcf, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "SNP VCF, duplicate keys (not bi-allelic)"

barcodes_list = []
tot_mtx_list = []
ad_mtx_list = []
for idx, dataset_id in enumerate(dataset_ids):
    barcodes = pd.read_table(sample_tsvs[idx], sep="\t", header=None, names=["BARCODE"])
    barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{dataset_id}"
    barcodes_list.append(barcodes)
    # cellsnp-lite emits one VCF + two MatrixMarket files per replicate
    tot_canon, ad_canon = map_allele_mat_to_snps(
        parent_keys,
        read_VCF(vcf_files[idx], addkey=True),
        mmread(tot_mtx_files[idx]).tocsr(),
        mmread(ad_mtx_files[idx]).tocsr(),
        len(barcodes),
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

all_barcodes = pd.concat(barcodes_list, axis=0, ignore_index=True)
cell_dataset_ids = np.repeat(
    np.arange(len(dataset_ids), dtype=np.int64),
    [len(b) for b in barcodes_list],
)
barcodes_full = pd.DataFrame(
    {
        "REP_ID": np.array(dataset_ids, dtype=str)[cell_dataset_ids],
        "BARCODE": all_barcodes["BARCODE"].to_numpy(),
    }
)
tot_mtx, ref_mtx, alt_mtx = hstack_replicate_mats(tot_mtx_list, ad_mtx_list)
a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, snps["PHASE"].to_numpy())

##################################################
num_snps_before = len(snps)
regions = read_BED(region_bed)

# feature_id (;-joined GTF genes) + feature_type, uniform across all assays
snps = annotate_feature_type(snps, gtf_file)

masks = [
    get_mask_by_region(snps, regions),
    get_mask_by_blacklist(snps, blacklist_bed),
]
if exon_only:
    masks.append(get_mask_by_exon(snps))
if is_rna_assay:
    # coverage filter only: RNA reads cover expressed genes, so drop SNPs outside
    # the h5ad feature set (feature_id itself stays GTF-derived from above)
    adata: sc.AnnData = sc.read_h5ad(h5ad_file)
    cov_mask = overlaps_any_range(snps, adata.var)
    logging.info(
        f"{assay_type} feature overlap: {np.sum(cov_mask)}/{len(snps)} "
        f"({np.sum(cov_mask) / len(snps):.3%})"
    )
    masks.append(cov_mask)

snps, snp_mask = apply_masks_to_df(snps, *masks)
snps = interp_pos_ranges(snps, regions, colname="region_id")

logging.info(f"#SNPs={np.sum(snp_mask)}/{num_snps_before} after filtering")

tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
a_mtx = a_mtx[snp_mask, :]
b_mtx = b_mtx[snp_mask, :]

sample_labels = [
    f"{dataset_id} {assay_type} {sample_type[0].upper()}"
    for dataset_id, sample_type in zip(dataset_ids, sample_types)
]

with PdfPages(out_qc_pdf) as pdf:
    plot_snp_depth(
        tot_mtx,
        sample_labels,
        qc_dir,
        f"{assay_type}.{run_id}",
        ref_mtx=ref_mtx,
        b_mtx=b_mtx,
        is_bulk=False,
        cell_dataset_ids=cell_dataset_ids,
        name_prefix="phase_and_concat",
        pdf=pdf,
        sample_id=sample_id,
    )
    plot_allele_freqs(
        snps,
        sample_labels,
        tot_mtx,
        ref_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        allele="ref",
        feature_label="SNP",
        suffix=".unphased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        sample_id=sample_id,
        pdf=pdf,
        cell_dataset_ids=cell_dataset_ids,
    )
    plot_allele_freqs(
        snps,
        sample_labels,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        allele="B",
        feature_label="SNP",
        suffix=".phased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        sample_id=sample_id,
        pdf=pdf,
        cell_dataset_ids=cell_dataset_ids,
    )

##################################################
logging.info("saving output files")
write_snp_info(snps, out_snp_info)
save_npz(out_tot_mtx_snp, tot_mtx)
save_npz(out_a_mtx_snp, a_mtx)
save_npz(out_b_mtx_snp, b_mtx)
snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
np.save(out_unique_snp_ids, snp_ids.to_numpy())
all_barcodes.to_csv(out_all_barcodes, sep="\t", header=False, index=False)
barcodes_full.to_csv(out_barcodes_full, sep="\t", header=True, index=False)
write_sample_ids(
    sample_id,
    dataset_ids,
    sample_types,
    [assay_type] * len(dataset_ids),
    out_sample_file,
)
logging.info("finished.")
