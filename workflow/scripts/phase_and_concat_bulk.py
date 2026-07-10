"""Joint phase-and-concat for ALL bulk assays on one shared SNP grid.

Every bulk replicate (across every bulk assay) is piled up against the same phased
het-SNP VCF, so ``canon_mat_one_replicate`` aligns them all to ONE shared parent SNP
set. We therefore build a single dense matrix directly (one pseudobulk column per
replicate, ordered assay-by-assay) instead of per-assay matrices that combine_counts
would have to re-union.

Parent het SNPs are kept when they are: in a region, not blacklisted, het-balanced in
EVERY normal pileup, and >= min_depth in EVERY sample. Outputs (under one per-stream
allele_dir/{bulkWGS,bulkWES}/ subdir) feed combine_counts directly.
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

from utils import *
from io_utils import *
from combine_counts_utils import *
from count_reads_utils import *
from aggregation_utils import *
from matplotlib.backends.backend_pdf import PdfPages
from plot_utils import plot_allele_freqs, plot_snp_depth_histogram


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


##################################################
log_file = snakemake_handle.log[0]
setup_logging(log_file)
logging.info("joint phase and concat for bulk assays")

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

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_name = snakemake_handle.params["sample_name"]
col_assays = list(snakemake_handle.params["col_assays"])
col_reps = list(snakemake_handle.params["col_reps"])
col_sample_types = list(snakemake_handle.params["col_sample_types"])
col_base_reps = list(snakemake_handle.params["col_base_reps"])
min_depth = int(snakemake_handle.params["min_depth"])
gamma = float(snakemake_handle.params["gamma"])
exon_only = snakemake_handle.params["exon_only"]
run_id = snakemake_handle.params["run_id"]

# outputs
snp_info = snakemake_handle.output["snp_info"]
out_tot_mtx_snp = snakemake_handle.output["tot_mtx_snp"]
out_a_mtx_snp = snakemake_handle.output["a_mtx_snp"]
out_b_mtx_snp = snakemake_handle.output["b_mtx_snp"]
out_sample_file = snakemake_handle.output["sample_file"]

n_samples = len(col_reps)
normal_cols = [k for k, st in enumerate(col_sample_types) if st == "normal"]
logging.info(
    f"sample_name={sample_name}, {n_samples} bulk samples across assays={col_assays}, "
    f"normal columns={normal_cols}"
)

##################################################
snps = read_VCF(snp_vcf, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "invalid bi-allelic SNP VCF file"

tot_mtx_list = []
ad_mtx_list = []
for idx in range(n_samples):
    tot_canon, ad_canon = canon_mat_one_replicate(
        parent_keys, vcf_files[idx], tot_mtx_files[idx], ad_mtx_files[idx], 1
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

tot_mtx, ref_mtx, alt_mtx = merge_mats(tot_mtx_list, ad_mtx_list)
a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, snps["PHASE"].to_numpy())

tot_mtx = tot_mtx.toarray()
ref_mtx = ref_mtx.toarray()
alt_mtx = alt_mtx.toarray()
a_mtx = a_mtx.toarray()
b_mtx = b_mtx.toarray()

for nc in normal_cols:
    log_ref_mapping_bias(
        ref_mtx[:, nc].astype(float),
        alt_mtx[:, nc].astype(float),
        label=f"Normal[{col_assays[nc]}:{col_reps[nc]}]",
    )

##################################################
num_snps_before = len(snps)
snp_mask = np.ones(len(snps), dtype=bool)
snp_mask, regions = apply_region_blacklist_masks(
    snps, snp_mask, region_bed, blacklist_bed
)

snps, _, _ = annotate_feature_type(snps, gtf_file)
snps.drop(columns=["gene_idx"], inplace=True, errors="ignore")

snp_mask &= get_mask_by_depth(snps, tot_mtx, min_dp=max(min_depth, 1))
for nc in normal_cols:
    snp_mask &= get_mask_by_het_balanced(snps, ref_mtx, alt_mtx, gamma, normal_idx=nc)

snp_mask = apply_exon_only_mask(snps, snp_mask, exon_only)

snps = snps.loc[snp_mask, :].reset_index(drop=True)
snps = assign_snp_bounderies(snps, regions, colname="region_id")

logging.info(f"#SNPs={np.sum(snp_mask)}/{num_snps_before} after filtering")

tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
a_mtx = a_mtx[snp_mask, :]
b_mtx = b_mtx[snp_mask, :]

##################################################
sample_labels = [f"{col_assays[k]}:{col_reps[k]}" for k in range(n_samples)]
plot_snp_depth_histogram(
    tot_mtx,
    sample_labels,
    qc_dir,
    f"bulk.{run_id}",
    ref_mtx=ref_mtx,
    is_bulk=True,
    cell_rep_idx=None,
    name_prefix="phase_and_concat",
)

af_pdf_path = snakemake_handle.output["qc_pdf"]
with PdfPages(af_pdf_path) as pdf:
    plot_allele_freqs(
        snps,
        sample_labels,
        tot_mtx,
        ref_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=False,
        allele="ref",
        unit="SNP",
        suffix=".unphased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        pdf=pdf,
    )
    plot_allele_freqs(
        snps,
        sample_labels,
        tot_mtx,
        b_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=False,
        allele="B",
        unit="SNP",
        suffix=".phased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        pdf=pdf,
    )

##################################################
logging.info("saving joint bulk output files")
out_cols = ["#CHR", "POS", "POS0", "START", "END", "GT", "PHASE"]
if "PS" in snps.columns:
    out_cols.append("PS")
out_cols += ["region_id", "feature_id", "feature_type"]
snps[out_cols].to_csv(snp_info, sep="\t", header=True, index=False)

np.savez_compressed(out_tot_mtx_snp, mat=tot_mtx)
np.savez_compressed(out_a_mtx_snp, mat=a_mtx)
np.savez_compressed(out_b_mtx_snp, mat=b_mtx)

sample_df = pd.DataFrame(
    {
        "SAMPLE": [f"{sample_name}_{rep}" for rep in col_reps],
        "SAMPLE_NAME": sample_name,
        "REP_ID": col_reps,
        "sample_type": col_sample_types,
        "assay_type": col_assays,
        "RDR_BASE_REP_ID": col_base_reps,
    }
)
sample_df.to_csv(out_sample_file, sep="\t", header=True, index=False)
logging.info("finished joint bulk phase_and_concat.")
