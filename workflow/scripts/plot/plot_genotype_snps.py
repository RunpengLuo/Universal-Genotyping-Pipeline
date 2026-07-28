"""Genotyping QC plot: reference-AF and depth histograms by genotype (het vs hom-alt)."""

import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def plot_genotype_af_depth(ref_af, total_depth, is_het, pdf, *, dpi=150):
    """One page: ref-AF and total-depth histograms, het vs hom-alt overlaid.

    ref_af : per-SNP reference allele fraction (ref / (ref + alt)).
    total_depth : per-SNP ref + alt allele depth.
    is_het : boolean mask, True for heterozygous SNPs (else homozygous ALT).
    """
    het = np.asarray(is_het, dtype=bool)
    ref_af = np.asarray(ref_af, dtype=float)
    total_depth = np.asarray(total_depth, dtype=float)

    fig, (ax_af, ax_dp) = plt.subplots(1, 2, figsize=(11, 4))
    groups = [
        ("het", "tab:blue", het),
        ("hom-alt", "tab:orange", ~het),
    ]
    for label, color, sel in groups:
        af = ref_af[sel & np.isfinite(ref_af)]
        if len(af):
            ax_af.hist(
                af,
                bins=50,
                range=(0, 1),
                alpha=0.6,
                color=color,
                label=f"{label} ({sel.sum()})",
            )
        dp = total_depth[sel & np.isfinite(total_depth)]
        if len(dp):
            hi = np.quantile(dp, 0.99) if len(dp) > 1 else dp.max()
            ax_dp.hist(
                dp[dp <= max(hi, 1)], bins=50, alpha=0.6, color=color, label=label
            )
    ax_af.axvline(0.5, color="grey", linestyle=":", linewidth=1)
    ax_af.set_xlim(0, 1)
    ax_af.set_xlabel("REF-allele frequency")
    ax_af.set_ylabel("# SNPs")
    ax_af.set_title(
        "REF-allele frequency by genotype (het vs hom-alt)",
        fontsize=9,
        fontweight="bold",
    )
    ax_af.legend(fontsize=8)
    ax_dp.set_xlabel("Total allele depth (ref + alt)")
    ax_dp.set_ylabel("# SNPs")
    ax_dp.set_title("Depth by genotype", fontsize=9, fontweight="bold")
    ax_dp.legend(fontsize=8)
    fig.tight_layout()
    pdf.savefig(fig, dpi=dpi)
    plt.close(fig)
