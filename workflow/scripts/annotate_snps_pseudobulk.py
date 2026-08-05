import os
import logging
import subprocess

snakemake_handle = snakemake

t = int(getattr(snakemake_handle, "threads", 1))
os.environ["OMP_NUM_THREADS"] = str(t)
os.environ["OPENBLAS_NUM_THREADS"] = str(t)
os.environ["MKL_NUM_THREADS"] = str(t)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(t)
os.environ["NUMEXPR_NUM_THREADS"] = str(t)

import numpy as np
import pandas as pd

from io_utils import read_VCF, get_chr_sizes
from utils import strip_chr_prefix

##################################################
"""
Given cellsnp-lite results from multiple replicates.
After annotation, only unique bi-allelic Het or Hom-Alt SNPs are kept.
"""
log_file = snakemake_handle.log[0]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
)

# inputs
raw_snp_vcfs = list(snakemake_handle.input["raw_snp_vcfs"])
genome_size = snakemake_handle.input["genome_size"]

# parameters
chroms = snakemake_handle.params["chroms"]
input_nochr = snakemake_handle.params["input_nochr"]
filter_nz_OTH = snakemake_handle.params["filter_nz_OTH"]
filter_hom_ALT = snakemake_handle.params["filter_hom_ALT"]
min_het_reads = int(snakemake_handle.params["min_het_reads"])
min_hom_dp = int(snakemake_handle.params["min_hom_dp"])
min_vaf_thres = float(snakemake_handle.params["min_vaf_thres"])
modalities = list(snakemake_handle.params["modalities"])

# outputs
snp_vcfs = list(snakemake_handle.output["snp_vcfs"])

logging.info(
    f"start annotate_snps_pseudobulk, filter_nz_OTH={filter_nz_OTH}, filter_hom_ALT={filter_hom_ALT}"
)
logging.info(
    f"min_het_reads={min_het_reads}, min_hom_dp={min_hom_dp}, min_vaf_thres={min_vaf_thres}"
)


def get_genotype(row):
    """Return the VCF GT string for a SNP row based on genotype boolean columns.

    Parameters
    ----------
    row : pd.Series
        Row with ``is_het``, ``is_hom_alt``, ``is_hom_ref`` boolean fields.

    Returns
    -------
    str
        One of ``"0/1"``, ``"1/1"``, ``"0/0"``, or ``"./."``.
    """
    if row["is_het"]:
        return "0/1"
    if row["is_hom_alt"]:
        return "1/1"
    if row["is_hom_ref"]:
        return "0/0"
    return "./."


KEY = ["#CHROM", "POS", "REF", "ALT"]
CNT = ["DP", "AD", "OTH"]
raw_snps_list = []
for idx, modality in enumerate(modalities):
    raw_snps = read_VCF(raw_snp_vcfs[idx], addkey=True)
    assert all(c in raw_snps.columns for c in KEY + CNT), "invalid cellsnp-lite format"
    for cnt in CNT:
        raw_snps[f"{cnt}{idx}"] = raw_snps[cnt].astype(np.int64)

    raw_snps = raw_snps[raw_snps["#CHROM"].astype(str).isin(chroms)]

    n_dup_rows = raw_snps.duplicated(subset="KEY", keep=False).sum()
    if n_dup_rows > 0:
        n_dup_keys = raw_snps.loc[
            raw_snps.duplicated(subset="KEY", keep=False), "KEY"
        ].nunique()
        logging.warning(
            f"{modality} have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
        )
        logging.warning("drop duplicated rows.")
        raw_snps = raw_snps.drop_duplicates(subset="KEY", keep="first").reset_index(
            drop=True
        )
    raw_snps_list.append(raw_snps)

base_snps = pd.concat(
    [raw_snps[["KEY", "#CHROM", "POS", "REF", "ALT"]] for raw_snps in raw_snps_list],
    ignore_index=True,
).drop_duplicates(subset=["#CHROM", "POS", "REF", "ALT"], keep="first")

dup_mask = base_snps.duplicated(subset="KEY", keep=False)
n_dup_rows = int(dup_mask.sum())
if n_dup_rows > 0:
    n_dup_keys = int(base_snps.loc[dup_mask, "KEY"].drop_duplicates().shape[0])
    logging.warning(
        f"merged df have {n_dup_rows} rows with duplicated positions={n_dup_keys}"
    )
    logging.warning("drop duplicated rows.")
    base_snps = base_snps.loc[~dup_mask, :]
base_snps = base_snps.reset_index(drop=True)
nsnps_before_genotyping = len(base_snps)

base_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort")
for cnt in CNT:
    base_snps[cnt] = 0
for idx, raw_snps in enumerate(raw_snps_list):
    right_cols = ["KEY"] + [f"{cnt}{idx}" for cnt in CNT]
    base_snps = pd.merge(
        left=base_snps,
        right=raw_snps[right_cols],
        on="KEY",
        how="left",
        validate="one_to_one",
        sort=False,
    )
    for cnt in CNT:
        base_snps[cnt] += base_snps[f"{cnt}{idx}"].fillna(0).astype(np.int64)
    base_snps = base_snps.drop(columns=[f"{cnt}{idx}" for cnt in CNT])

base_snps["REF_COUNT"] = base_snps["DP"] - base_snps["AD"]
base_snps["AR"] = base_snps["AD"] / base_snps["DP"]

base_snps["is_het"] = (
    base_snps[["AD", "REF_COUNT"]].min(axis=1) >= min_het_reads
) & base_snps["AR"].between(min_vaf_thres, 1 - min_vaf_thres)

base_snps["is_hom_alt"] = (base_snps["AD"] == base_snps["DP"]) & (
    base_snps["DP"] >= min_hom_dp
)
base_snps["is_hom_ref"] = (base_snps["AD"] == 0) & (base_snps["DP"] >= min_hom_dp)

base_snps["SAMPLE"] = base_snps.apply(get_genotype, axis=1)

chrom_sizes = {
    f"chr{strip_chr_prefix(k)}": v for k, v in get_chr_sizes(genome_size).items()
}

keep_gts = ["0/1"]
if not filter_hom_ALT:
    keep_gts.append("1/1")
base_snps = base_snps[base_snps["SAMPLE"].isin(keep_gts)].reset_index(drop=True)

nz_oth = base_snps["OTH"] > 0
nnz_oth = np.sum(nz_oth)
logging.info(f"#nz-OTH SNPs={nnz_oth}/{nsnps_before_genotyping}")
if filter_nz_OTH:
    logging.info("SNPs with nonzero OTHs are filtered.")
    base_snps = base_snps[~nz_oth]

num_kept = len(base_snps)
logging.info(f"#kept SNPs={num_kept}/{nsnps_before_genotyping}")

final_snps = base_snps.sort_values(["#CHROM", "POS"], kind="mergesort").reset_index(
    drop=True
)

final_snps["FORMAT"] = "GT"
final_snps["ID"] = "."
final_snps["QUAL"] = "."
final_snps["FILTER"] = "PASS"
final_snps["INFO"] = (
    "AD="
    + final_snps["AD"].astype(str)
    + ";DP="
    + final_snps["DP"].astype(str)
    + ";OTH="
    + final_snps["OTH"].astype(str)
)

cols = [
    "#CHROM",
    "POS",
    "ID",
    "REF",
    "ALT",
    "QUAL",
    "FILTER",
    "INFO",
    "FORMAT",
    "SAMPLE",
]
final_snps = final_snps[cols]

final_snps_chs = final_snps.groupby("#CHROM", sort=False)
for chrom, out_snp_vcf in zip(chroms, snp_vcfs):
    contig = strip_chr_prefix(chrom) if input_nochr else chrom
    chrom_length = chrom_sizes[chrom]
    with open(out_snp_vcf[:-3], "w") as fd:
        fd.writelines(
            [
                "##fileformat=VCFv4.2\n",
                f"##contig=<ID={contig},length={chrom_length}>\n",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Pseudobulk genotype">\n',
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth, REF+ALT">\n',
                '##INFO=<ID=AD,Number=1,Type=Integer,Description="Allele Depth for ALT allele">\n',
                '##INFO=<ID=OTH,Number=1,Type=Integer,Description="Allele Depth other than REF and ALT">\n',
            ]
        )
        fd.write("\t".join(cols) + "\n")
        if chrom in final_snps_chs.groups:
            out = final_snps_chs.get_group(chrom)[cols].copy()
            out["#CHROM"] = contig
            out.to_csv(fd, sep="\t", index=False, header=False)

    # out_snp_vcf[:-3] will be removed by bgzip.
    subprocess.run(["bgzip", "-f", out_snp_vcf[:-3]], check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out_snp_vcf], check=True)
logging.info("finished annotate_snps_pseudobulk")
