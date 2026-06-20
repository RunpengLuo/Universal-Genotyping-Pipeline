"""Phase-and-concat allele counts for a single non-bulk assay (per assay_type).

Merges the assay's replicate pileups onto the shared parent SNP set, phases ref/alt
into A/B, filters SNPs (region, blacklist, RNA gene overlap / ATAC none), and writes
sparse (SNPs x cells) matrices plus per-cell barcode bookkeeping. Joint binning across
the sample's assays happens downstream in combine_counts_nonbulk.
"""

import os
import logging

snakemake_handle = snakemake

t = int(getattr(snakemake_handle, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd
from scipy.sparse import save_npz

from utils import *
from io_utils import *
from combine_counts_utils import *
from count_reads_utils import *
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_allele_freqs, plot_snp_depth_histogram
from aggregation_utils import *


def annotate_feature_type(snps, gtf_file):
    """Annotate SNPs with feature_type (exon/intron/intergenic) using a GTF file."""
    genes_gtf = read_genes_gtf_file(gtf_file, id_col="gene_id")[
        ["gene_id", "#CHR", "START", "END"]
    ]
    genes_gtf["gene_idx"] = np.arange(len(genes_gtf))
    snps = assign_pos_to_range(snps, genes_gtf, ref_id="gene_idx", pos_col="POS0")
    gene_mask = snps["gene_idx"].notna()

    exons_gtf = read_exons_gtf_file(gtf_file)
    exons_gtf["exon_idx"] = np.arange(len(exons_gtf))
    snps = assign_pos_to_range(snps, exons_gtf, ref_id="exon_idx", pos_col="POS0")

    snps["feature_type"] = "intergenic"
    snps.loc[gene_mask, "feature_type"] = "intron"
    snps.loc[snps["exon_idx"].notna(), "feature_type"] = "exon"

    snps.drop(columns=["exon_idx"], inplace=True, errors="ignore")
    return snps, genes_gtf, gene_mask


##################################################
log_file = snakemake_handle.log[0]
setup_logging(log_file)
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
sample_name = snakemake_handle.params["sample_name"]
assay_type = snakemake_handle.params["assay_type"]
rep_ids = snakemake_handle.params["rep_ids"]
sample_types = snakemake_handle.params["sample_types"]
exon_only = snakemake_handle.params["exon_only"]
run_id = snakemake_handle.params["run_id"]

# outputs
snp_info = snakemake_handle.output["snp_info"]
tot_mtx_snp = snakemake_handle.output["tot_mtx_snp"]
a_mtx_snp = snakemake_handle.output["a_mtx_snp"]
b_mtx_snp = snakemake_handle.output["b_mtx_snp"]
unique_snp_ids = snakemake_handle.output["unique_snp_ids"]
cell_snp_Aallele = snakemake_handle.output["cell_snp_Aallele"]
cell_snp_Ballele = snakemake_handle.output["cell_snp_Ballele"]
out_all_barcodes = snakemake_handle.output["all_barcodes"]
out_barcodes_full = snakemake_handle.output["barcodes_full"]
sample_file = snakemake_handle.output["sample_file"]

is_rna_assay = ASSAY_TYPE2MODALITY[assay_type] == "RNA"

##################################################
logging.info(f"sample_name={sample_name}, assay_type={assay_type}, rep_ids={rep_ids}")

snps = read_VCF(snp_vcf, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

barcodes_list = []
tot_mtx_list = []
ad_mtx_list = []
for idx, rep_id in enumerate(rep_ids):
    barcodes = pd.read_table(sample_tsvs[idx], sep="\t", header=None, names=["BARCODE"])
    barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{rep_id}"
    barcodes_list.append(barcodes)
    tot_canon, ad_canon = canon_mat_one_replicate(
        parent_keys, vcf_files[idx], tot_mtx_files[idx], ad_mtx_files[idx], len(barcodes)
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

all_barcodes = pd.concat(barcodes_list, axis=0, ignore_index=True)
cell_rep_idx = np.repeat(
    np.arange(len(rep_ids), dtype=np.int64),
    [len(b) for b in barcodes_list],
)
barcodes_full = pd.DataFrame(
    {
        "REP_ID": np.array(rep_ids, dtype=str)[cell_rep_idx],
        "BARCODE": all_barcodes["BARCODE"].to_numpy(),
    }
)
tot_mtx, ref_mtx, alt_mtx = merge_mats(tot_mtx_list, ad_mtx_list)
a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, snps["PHASE"].to_numpy())

##################################################
num_snps_before = len(snps)

snp_mask = np.ones(len(snps), dtype=bool)
regions = read_region_file(region_bed)
region_mask = get_mask_by_region(snps, regions)
logging.info(f"region filter: {np.sum(region_mask)}/{len(snps)} SNPs passed")
snp_mask &= region_mask

if blacklist_bed is not None:
    bl_regions = read_region_file(blacklist_bed)
    bl_mask = get_mask_by_region(snps, bl_regions)
    logging.info(f"blacklist filter: {np.sum(bl_mask)}/{len(snps)} SNPs in blacklist")
    snp_mask &= ~bl_mask

snps, genes_gtf, gene_mask = annotate_feature_type(snps, gtf_file)
snps.drop(columns=["gene_idx"], inplace=True, errors="ignore")

if is_rna_assay:
    # filter SNPs to those overlapping an RNA feature (gene) from the h5ad
    logging.info("annotate SNPs with feature_id")
    adata: sc.AnnData = sc.read_h5ad(h5ad_file)
    feature_df = adata.var.reset_index(drop=False).rename(
        columns={"index": "feature_id"}
    )
    feature_df["feature_idx"] = np.arange(len(feature_df))

    snps = assign_pos_to_range(snps, feature_df, ref_id="feature_idx", pos_col="POS0")
    snp_mask &= snps["feature_idx"].notna().to_numpy()
    logging.info(
        f"{assay_type} feature overlap: {np.sum(snp_mask)}/{len(snps)} "
        f"({np.sum(snp_mask) / len(snps):.3%})"
    )
else:
    # scATAC: no feature-level filtering; region filter (above) already applied
    snps["feature_id"] = "intergenic"

_n_exon = int((snps["feature_type"] == "exon").sum())
logging.info(f"#exonic SNPs: {_n_exon}/{len(snps)} ({_n_exon / max(len(snps), 1):.3%})")
if exon_only:
    exon_mask = (snps["feature_type"] == "exon").to_numpy()
    logging.info(f"exon filter: {np.sum(exon_mask)}/{len(snps)} SNPs passed")
    snp_mask &= exon_mask

snps = snps.loc[snp_mask, :].reset_index(drop=True)
if is_rna_assay:
    snps["feature_idx"] = snps["feature_idx"].astype(feature_df["feature_idx"].dtype)
    snps = pd.merge(
        left=snps,
        right=feature_df[["feature_idx", "feature_id"]],
        how="left",
        on="feature_idx",
        sort=False,
    ).reset_index(drop=True)
# TODO refine by gene interval?
snps["START"] = snps["POS0"]
snps["END"] = snps["POS"]

snps = assign_snp_bounderies(snps, regions, colname="region_id")

logging.info(f"#SNPs={np.sum(snp_mask)}/{num_snps_before} after filtering")

tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
a_mtx = a_mtx[snp_mask, :]
b_mtx = b_mtx[snp_mask, :]

plot_snp_depth_histogram(
    tot_mtx,
    rep_ids,
    qc_dir,
    f"{assay_type}.{run_id}",
    ref_mtx=ref_mtx,
    is_bulk=False,
    cell_rep_idx=cell_rep_idx,
    name_prefix="phase_and_concat",
)

af_pdf_path = os.path.join(qc_dir, f"phase_and_concat.snp_allele_freq.{assay_type}.{run_id}.pdf")
with PdfPages(af_pdf_path) as pdf:
    plot_allele_freqs(
        snps,
        rep_ids,
        tot_mtx,
        ref_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        allele="ref",
        unit="SNP",
        suffix=".unphased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        pdf=pdf,
        cell_rep_idx=cell_rep_idx,
    )
    plot_allele_freqs(
        snps,
        rep_ids,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=True,
        allele="B",
        unit="SNP",
        suffix=".phased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        pdf=pdf,
        cell_rep_idx=cell_rep_idx,
    )

##################################################
logging.info("saving output files")
snps[
    [
        "#CHR",
        "POS",
        "POS0",
        "START",
        "END",
        "GT",
        "PHASE",
        "region_id",
        "feature_id",
        "feature_type",
    ]
].to_csv(snp_info, sep="\t", header=True, index=False)
save_npz(tot_mtx_snp, tot_mtx)
save_npz(a_mtx_snp, a_mtx)
save_npz(b_mtx_snp, b_mtx)
snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
np.save(unique_snp_ids, snp_ids.to_numpy())
save_npz(cell_snp_Aallele, a_mtx)
save_npz(cell_snp_Ballele, b_mtx)
all_barcodes.to_csv(out_all_barcodes, sep="\t", header=False, index=False)
barcodes_full.to_csv(out_barcodes_full, sep="\t", header=True, index=False)
sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df["sample_type"] = sample_types
sample_df.to_csv(sample_file, sep="\t", header=True, index=False)
logging.info("finished.")
