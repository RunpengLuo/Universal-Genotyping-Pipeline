"""Genotyping QC: het vs hom-alt allele-fraction diagnostic (bulk genotyped VCFs).

Runpeng Luo (2026-07-27)

Reads the per-chromosome genotyped VCFs (``snps/chr{c}.vcf.gz``), which retain het +
hom-alt SNPs with the genotyping sample's ``FORMAT/AD`` and ``DP``, and renders a QC PDF:
a genome-wide reference allele fraction (AF) scatter colored by genotype (het vs hom-alt),
plus per-genotype AF and depth histograms. hom-ref is never genotyped, so only het vs
hom-alt are shown.

Inputs (snakemake.input):
    vcfs: per-chromosome genotyped VCFs (het + hom-alt, with FORMAT/AD, DP).
    genome_size: chrom sizes file for genome-wide coordinates.
Outputs (snakemake.output):
    qc_pdf: genotype SNP QC PDF.
"""

import logging

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import read_VCF
from plot_genome import plot_1d_sample
from plot_genotype_snps import plot_genotype_af_depth
from utils import setup_logging

_HET_GTS = {"0/1", "1/0", "0|1", "1|0"}


def _parse_ad(ad_series):
    """Parse a FORMAT/AD string column ('ref,alt') into ref and alt float arrays."""
    parts = ad_series.fillna("").str.split(",", n=2, expand=True)
    ref = pd.to_numeric(parts[0], errors="coerce").to_numpy(dtype=float)
    if parts.shape[1] > 1:
        alt = pd.to_numeric(parts[1], errors="coerce").to_numpy(dtype=float)
    else:
        alt = np.full(len(ref), np.nan)
    return ref, alt


setup_logging(snakemake.log[0])
vcfs = snakemake.input.vcfs
genome_size = snakemake.input.genome_size
out_pdf = snakemake.output.qc_pdf

frames = [df for vcf in vcfs if (df := read_VCF(vcf, addkey=True)) is not None]
if not frames:
    logging.warning("no genotyped SNPs found; writing placeholder QC PDF")
    with PdfPages(out_pdf) as pdf:
        fig = plt.figure(figsize=(8, 2))
        fig.text(0.5, 0.5, "no genotyped SNPs", ha="center", va="center")
        pdf.savefig(fig)
        plt.close(fig)
else:
    snps = pd.concat(frames, ignore_index=True)
    if "AD" not in snps.columns:
        raise ValueError(
            "genotyped VCF has no FORMAT/AD; cannot compute allele fraction"
        )

    ref, alt = _parse_ad(snps["AD"])
    total = ref + alt
    with np.errstate(invalid="ignore", divide="ignore"):
        ref_af = np.where(total > 0, ref / total, np.nan)
    is_het = snps["GT"].isin(_HET_GTS).to_numpy()

    logging.info(
        f"genotype QC: {len(snps)} SNPs = {int(is_het.sum())} het + "
        f"{int((~is_het).sum())} hom-alt; median ref-AF "
        f"het={np.nanmedian(ref_af[is_het]):.3f}, "
        f"hom-alt={np.nanmedian(ref_af[~is_het]):.3f}"
    )

    pos_df = snps[["#CHR", "POS"]].copy()
    pos_df["#CHR"] = pos_df["#CHR"].astype(str)
    with PdfPages(out_pdf) as pdf:
        plot_1d_sample(
            pos_df,
            ref_af,
            genome_size,
            out_pdf,
            unit="SNP",
            val_type="AF",
            mask=is_het,
            mask_labels=("het", "hom-alt"),
            mask_colors=("tab:blue", "tab:orange"),
            pdf=pdf,
        )
        plot_genotype_af_depth(ref_af, total, is_het, pdf)
    logging.info(f"saved genotype SNP QC to {out_pdf}")
