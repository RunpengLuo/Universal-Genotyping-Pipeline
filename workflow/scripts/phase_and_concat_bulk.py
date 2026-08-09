"""Joint phase-and-concat for ALL bulk assays on one shared SNP set.

Every bulk replicate (across every bulk assay) is piled up against the same phased
het-SNP VCF, so ``map_allele_mat_to_snps`` aligns them all to ONE shared parent SNP
set. We therefore build a single dense matrix directly (one pseudobulk column per
replicate, ordered assay-by-assay) instead of per-assay matrices that combine_counts
would have to re-union.

Parent het SNPs are kept when they are: in a region, not blacklisted, het-balanced in
EVERY normal pileup, and >= min_depth in EVERY sample. Outputs (under one per-stream
allele_dir/{bulkWGS,bulkWES}/ subdir) feed combine_counts directly.
"""

import logging

snakemake_handle = snakemake

from utils import set_omp_threads, setup_logging, log_hist, maybe_path

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd

from io_utils import (
    read_BED,
    read_VCF,
    read_bcftools_pileup_counts,
    write_sample_ids,
    write_snp_info,
)
from combine_counts_utils import (
    apply_masks_to_df,
    build_pos_ranges,
    get_mask_by_blacklist,
    get_mask_by_depth,
    get_mask_by_exon,
    get_mask_by_het_balanced,
    get_mask_by_region,
    hstack_replicate_mats,
    map_allele_mat_to_snps,
)
from phasing_utils import apply_phase_to_mat
from feature_utils import annotate_feature_type
from matplotlib.backends.backend_pdf import PdfPages
from plot_alleles import plot_allele_freqs, plot_snp_depth
from plot_utils import observation_order


##################################################
logging.info("joint phase and concat for bulk assays")

# inputs
counts_files = snakemake_handle.input["counts"]
snp_vcf = snakemake_handle.input["snp_vcf"]
region_bed = snakemake_handle.input["region_bed"]
genome_size = snakemake_handle.input["genome_size"]
gtf_file = maybe_path(snakemake_handle.input["gtf_file"])
blacklist_bed = maybe_path(snakemake_handle.input["blacklist_bed"])

# parameters
qc_dir = snakemake_handle.params["qc_dir"]
sample_id = snakemake_handle.params["sample_id"]
dataset_assays = list(snakemake_handle.params["dataset_assays"])
dataset_ids = list(snakemake_handle.params["dataset_ids"])
sample_types = list(snakemake_handle.params["sample_types"])
base_dataset_ids = list(snakemake_handle.params["base_dataset_ids"])
min_depth = int(snakemake_handle.params["min_depth"])
gamma = float(snakemake_handle.params["gamma"])
exon_only = snakemake_handle.params["exon_only"]
run_id = snakemake_handle.params["run_id"]

# outputs
out_snp_info = snakemake_handle.output["snp_info"]
out_tot_mtx_snp = snakemake_handle.output["tot_mtx_snp"]
out_a_mtx_snp = snakemake_handle.output["a_mtx_snp"]
out_b_mtx_snp = snakemake_handle.output["b_mtx_snp"]
out_sample_file = snakemake_handle.output["sample_file"]
out_qc_pdf = snakemake_handle.output["qc_pdf"]

n_samples = len(dataset_ids)
normal_obs = [k for k, st in enumerate(sample_types) if st == "normal"]
logging.info(
    f"sample_id={sample_id}, {n_samples} bulk samples across assays={dataset_assays}, "
    f"normal columns={normal_obs}"
)

##################################################
snps = read_VCF(snp_vcf, addkey=True, add_phase1=True, add_pos0=True)
parent_keys = pd.Index(snps["KEY"])
assert not parent_keys.duplicated().any(), "SNP VCF, duplicate keys (not bi-allelic)"
parent_alt_by_key = dict(zip(snps["KEY"], snps["ALT"]))

tot_mtx_list = []
ad_mtx_list = []
for idx in range(n_samples):
    child_snps, tot_child, ad_child = read_bcftools_pileup_counts(
        counts_files[idx], parent_alt_by_key
    )
    tot_canon, ad_canon = map_allele_mat_to_snps(
        parent_keys, child_snps, tot_child, ad_child, 1
    )
    tot_mtx_list.append(tot_canon)
    ad_mtx_list.append(ad_canon)

tot_mtx, ref_mtx, alt_mtx = hstack_replicate_mats(tot_mtx_list, ad_mtx_list)
a_mtx, b_mtx = apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, snps["PHASE"].to_numpy())

tot_mtx = tot_mtx.toarray()
ref_mtx = ref_mtx.toarray()
alt_mtx = alt_mtx.toarray()
a_mtx = a_mtx.toarray()
b_mtx = b_mtx.toarray()

# REF/(REF+ALT) per normal: mass sitting off 0.5 is reference mapping bias
for nc in normal_obs:
    total = (ref_mtx[:, nc] + alt_mtx[:, nc]).astype(float)
    covered = total > 0
    log_hist(
        ref_mtx[covered, nc] / total[covered],
        f"Normal[{dataset_assays[nc]}:{dataset_ids[nc]}] REF/(REF+ALT) over {len(total)} SNPs",
    )

##################################################
num_snps_before = len(snps)
regions = read_BED(region_bed)
snps = annotate_feature_type(snps, gtf_file)

masks = [
    get_mask_by_region(snps, regions),
    get_mask_by_blacklist(snps, blacklist_bed),
    get_mask_by_depth(snps, tot_mtx, min_dp=max(min_depth, 1)),
    *(
        get_mask_by_het_balanced(snps, ref_mtx, alt_mtx, gamma, normal_idx=nc)
        for nc in normal_obs
    ),
]
if exon_only:
    masks.append(get_mask_by_exon(snps))

snps, snp_mask = apply_masks_to_df(snps, *masks)
snps = build_pos_ranges(snps, regions, colname="region_id")

logging.info(f"#SNPs={np.sum(snp_mask)}/{num_snps_before} after filtering")

tot_mtx = tot_mtx[snp_mask, :]
ref_mtx = ref_mtx[snp_mask, :]
a_mtx = a_mtx[snp_mask, :]
b_mtx = b_mtx[snp_mask, :]

##################################################
sample_labels = [f"{dataset_assays[k]}:{dataset_ids[k]}" for k in range(n_samples)]
obs_order = observation_order(dataset_assays, sample_types, dataset_ids)

with PdfPages(out_qc_pdf) as pdf:
    plot_snp_depth(
        tot_mtx,
        sample_labels,
        qc_dir,
        f"bulk.{run_id}",
        ref_mtx=ref_mtx,
        b_mtx=b_mtx,
        is_bulk=True,
        cell_dataset_ids=None,
        name_prefix="phase_and_concat",
        pdf=pdf,
        obs_order=obs_order,
        sample_id=sample_id,
    )
    plot_allele_freqs(
        snps,
        sample_labels,
        tot_mtx,
        ref_mtx,
        genome_size,
        qc_dir,
        apply_pseudobulk=False,
        allele="ref",
        feature_label="SNP",
        suffix=".unphased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        sample_id=sample_id,
        obs_order=obs_order,
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
        feature_label="SNP",
        suffix=".phased",
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        run_id=run_id,
        sample_id=sample_id,
        obs_order=obs_order,
        pdf=pdf,
    )

##################################################
logging.info("saving phased allele count mats to files")
write_snp_info(snps, out_snp_info)

np.savez_compressed(out_tot_mtx_snp, mat=tot_mtx)
np.savez_compressed(out_a_mtx_snp, mat=a_mtx)
np.savez_compressed(out_b_mtx_snp, mat=b_mtx)

write_sample_ids(
    sample_id,
    dataset_ids,
    sample_types,
    dataset_assays,
    out_sample_file,
    rdr_base_dataset_ids=base_dataset_ids,
)
logging.info("finished joint bulk phase_and_concat.")
