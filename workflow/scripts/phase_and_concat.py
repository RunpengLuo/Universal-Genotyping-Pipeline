import os, logging

t = int(getattr(snakemake, "threads", 1))
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
from switchprobs import *


def log_ref_mapping_bias(ref_counts, alt_counts, label=""):
    """Log REF/(REF+ALT) summary stats to detect reference mapping bias."""
    total = ref_counts + alt_counts
    pos = total > 0
    n_pos = int(np.sum(pos))
    logging.info(f"{label}: {n_pos}/{len(total)} SNPs with total > 0")
    if n_pos > 0:
        ratio = ref_counts[pos] / total[pos]
        logging.info(
            f"{label} REF/(REF+ALT) stats: "
            f"min={np.min(ratio):.4f}, max={np.max(ratio):.4f}, "
            f"median={np.median(ratio):.4f}, mean={np.mean(ratio):.4f}"
        )


def annotate_feature_type(snps, gtf_file):
    """Annotate SNPs with feature_type (exon/intron/intergenic) using a GTF file.

    Returns the gene_mask (bool Series) indicating which SNPs fall within a gene.
    """
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
setup_logging(snakemake.log[0])
logging.info("phase and concat allele-level count matrices")

vcf_files = snakemake.input["vcfs"]
sample_tsvs = snakemake.input["sample_tsvs"]
tot_mtx_files = snakemake.input["tot_mtxs"]
ad_mtx_files = snakemake.input["ad_mtxs"]
snp_vcf = snakemake.input["snp_vcf"]
qc_dir = snakemake.params["qc_dir"]
qc_prefix = "phase_and_concat"
os.makedirs(qc_dir, exist_ok=True)
run_id = getattr(snakemake.params, "run_id", "")
qc_stamp = ".".join(p for p in (snakemake.params["assay_type"], run_id) if p)

region_bed = snakemake.input["region_bed"]
genome_size = snakemake.input["genome_size"]
gtf_file = maybe_path(snakemake.input["gtf_file"])
blacklist_bed = maybe_path(snakemake.input["blacklist_bed"])

sample_name = snakemake.params["sample_name"]
assay_type = snakemake.params["assay_type"]
rep_ids = snakemake.params["rep_ids"]
sample_types = snakemake.params["sample_types"]

is_bulk_assay = assay_type in BULK_ASSAYS

##################################################
logging.info(
    f"sample_name={sample_name}, assay_type={assay_type}, "
    f"is_bulk_assay={is_bulk_assay}, rep_ids={rep_ids}"
)

snps = read_VCF(snp_vcf, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

if is_bulk_assay:
    has_normal = "normal" in sample_types

barcodes_list = []
tot_mtx_list = []
ad_mtx_list = []
for idx, rep_id in enumerate(rep_ids):
    barcodes = pd.read_table(sample_tsvs[idx], sep="\t", header=None, names=["BARCODE"])
    barcodes["BARCODE"] = barcodes["BARCODE"].astype(str) + f"_{rep_id}"
    nbarcodes = len(barcodes)
    barcodes_list.append(barcodes)
    if is_bulk_assay:
        nbarcodes = 1
    tot_canon, ad_canon = canon_mat_one_replicate(
        parent_keys, vcf_files[idx], tot_mtx_files[idx], ad_mtx_files[idx], nbarcodes
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

all_barcodes = pd.concat(barcodes_list, axis=0, ignore_index=True)
if is_bulk_assay:
    cell_rep_idx = None
    barcodes_full = None
else:
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

if is_bulk_assay:
    tot_mtx = tot_mtx.toarray()
    ref_mtx = ref_mtx.toarray()
    alt_mtx = alt_mtx.toarray()
    a_mtx = a_mtx.toarray()
    b_mtx = b_mtx.toarray()

    if has_normal:
        log_ref_mapping_bias(
            ref_mtx[:, 0].astype(float),
            alt_mtx[:, 0].astype(float),
            label="Normal",
        )

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
if is_bulk_assay:
    snps["feature_id"] = "intergenic"
    if gene_mask.any():
        snps.loc[gene_mask, "gene_idx"] = snps.loc[gene_mask, "gene_idx"].astype(int)
        snps.loc[gene_mask, "feature_id"] = (
            genes_gtf.set_index("gene_idx")
            .loc[snps.loc[gene_mask, "gene_idx"], "gene_id"]
            .values
        )

    snps.drop(columns=["gene_idx"], inplace=True, errors="ignore")

    snp_mask &= get_mask_by_depth(
        snps, tot_mtx, min_dp=max(int(snakemake.params["min_depth"]), 1)
    )
    if has_normal:
        normal_mask = get_mask_by_het_balanced(
            snps, ref_mtx, alt_mtx, float(snakemake.params["gamma"]), normal_idx=0
        )
        snp_mask &= normal_mask
else:
    snps.drop(columns=["gene_idx"], inplace=True, errors="ignore")

    logging.info("annotate SNPs with feature_id")
    adata: sc.AnnData = sc.read_h5ad(snakemake.input["h5ad_file"])
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

_n_exon = int((snps["feature_type"] == "exon").sum())
_n_total = len(snps)
logging.info(f"#exonic SNPs: {_n_exon}/{_n_total} ({_n_exon / max(_n_total, 1):.3%})")

if snakemake.params["exon_only"]:
    exon_mask = (snps["feature_type"] == "exon").to_numpy()
    logging.info(f"exon filter: {np.sum(exon_mask)}/{len(snps)} SNPs passed")
    snp_mask &= exon_mask

snps = snps.loc[snp_mask, :].reset_index(drop=True)
if not is_bulk_assay:
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

num_snps_after = np.sum(snp_mask)
logging.info(f"#SNPs={num_snps_after}/{num_snps_before} after filtering")

tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
alt_mtx = alt_mtx[snp_mask, :]
a_mtx = a_mtx[snp_mask, :]
b_mtx = b_mtx[snp_mask, :]

plot_snp_depth_histogram(
    tot_mtx,
    rep_ids,
    qc_dir,
    qc_stamp,
    ref_mtx=ref_mtx,
    is_bulk=is_bulk_assay,
    cell_rep_idx=cell_rep_idx,
    name_prefix=qc_prefix,
)

af_pdf_path = qc_path(qc_dir, qc_prefix, "snp_allele_freq.pdf", qc_stamp)
with PdfPages(af_pdf_path) as pdf:
    plot_allele_freqs(
        snps,
        rep_ids,
        tot_mtx,
        ref_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=not is_bulk_assay,
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
        apply_pseudobulk=not is_bulk_assay,
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
].to_csv(snakemake.output["snp_info"], sep="\t", header=True, index=False)
if is_bulk_assay:
    np.savez_compressed(snakemake.output["tot_mtx_snp"], mat=tot_mtx)
    np.savez_compressed(snakemake.output["a_mtx_snp"], mat=a_mtx)
    np.savez_compressed(snakemake.output["b_mtx_snp"], mat=b_mtx)
else:
    save_npz(snakemake.output["tot_mtx_snp"], tot_mtx)
    save_npz(snakemake.output["a_mtx_snp"], a_mtx)
    save_npz(snakemake.output["b_mtx_snp"], b_mtx)
    snp_ids = snps["#CHR"].astype(str) + "_" + snps["POS"].astype(str)
    np.save(snakemake.output["unique_snp_ids"], snp_ids.to_numpy())
    save_npz(snakemake.output["cell_snp_Aallele"], a_mtx)
    save_npz(snakemake.output["cell_snp_Ballele"], b_mtx)
all_barcodes.to_csv(snakemake.output["all_barcodes"], sep="\t", header=False, index=False)
if barcodes_full is not None:
    barcodes_full.to_csv(
        snakemake.output["barcodes_full"], sep="\t", header=True, index=False
    )
sample_df = pd.DataFrame({"SAMPLE": [f"{sample_name}_{rep_id}" for rep_id in rep_ids]})
sample_df["SAMPLE_NAME"] = sample_name
sample_df["REP_ID"] = rep_ids
sample_df["sample_type"] = sample_types
sample_df.to_csv(snakemake.output["sample_file"], sep="\t", header=True, index=False)
logging.info("finished.")
