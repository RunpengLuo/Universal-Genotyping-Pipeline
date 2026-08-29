"""Shared plotting base: genome axis, shading, page layout, histograms.

Last update: 2026-08-11

Functions:
- _get_axis, _shade, _load_shading: the cnplot genome axis and region shading
- _observation_labels: the one place a per-observation display label is built
- _suptitle, _finish_page, _bold_chrnames: page titling and layout
- _hist_with_stats: a histogram with a mean/median-annotated title
- _val_full: full title name for a value-type abbreviation
"""

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cnplot import GenomeAxis, read_bed, shade_regions

from io_utils import read_chrom_sizes


# Full names for value-type abbreviations, used in figure titles only (ylabels keep
# the short form). AF is the reference-allele frequency (ref / total).
_VAL_TYPE_FULL = {
    "AF": "REF-allele frequency",
    "BAF": "B-allele frequency",
    "RD": "Read depth",
    "RDR": "Read-depth ratio",
}


def _observation_labels(dataset_ids, assay_types, sample_types):
    """One display label per observation: ``"{dataset_id} {assay_type} {T|N}"``.

    The single place this string is built. Plot functions take the three identifying
    columns and compose it here, so no caller assembles display text.

    Args:
        dataset_ids: Dataset id per observation.
        assay_types: Assay type per observation; pass ``[assay_type] * n`` for a
            single-assay frame.
        sample_types: ``tumor``/``normal`` per observation; only its first letter is
            shown, uppercased.

    Returns:
        List of labels, one per observation, in the given order.
    """
    return [
        f"{dataset_id} {assay_type} {str(sample_type)[0].upper()}"
        for dataset_id, assay_type, sample_type in zip(
            dataset_ids, assay_types, sample_types
        )
    ]


def _val_full(val_type):
    """Full title name for a value-type abbreviation (falls back to the input)."""
    return _VAL_TYPE_FULL.get(val_type, val_type)


# ---------------------------------------------------------------------------
# cnplot genome axis + region shading (shared by the 1D scatter plots)
# ---------------------------------------------------------------------------

_AXIS_CACHE = {}


def _get_axis(genome_size, chroms):
    """Whole-genome cnplot GenomeAxis restricted to the chromosomes in the data.

    Full chromosome widths (``collapse_gaps=False``); region/blacklist overlays are
    a separate ``shade_regions`` concern. Cached by (``genome_size``, chromosome
    set) since the axis depends only on the reference.
    """
    chroms = list(dict.fromkeys(str(c) for c in chroms))
    key = (genome_size, frozenset(chroms))
    axis = _AXIS_CACHE.get(key)
    if axis is None:
        keep = set(chroms)
        excluded = [c for c in read_chrom_sizes(genome_size) if c not in keep]
        axis = GenomeAxis(
            None, genome_size, excluded_chroms=excluded, collapse_gaps=False
        )
        _AXIS_CACHE[key] = axis
    return axis


def _shade(ax, axis, region_df, blacklist_df):
    """Shade masked (blacklist) regions gray; callable regions keep the white background.

    *region_df* (the callable whitelist) is unshaded so the plotted background stays
    white; it is kept in the signature for caller stability. Blacklist intervals are
    unioned first so overlaps do not compound the alpha.
    """

    def _merge_ranges(df):
        """Union overlapping/touching [START, END) ranges per ``#CHR``.

        Shading each raw range separately stacks the alpha where they overlap, so
        overlapping masked regions render darker than a single range. Merging first
        keeps the fill uniform.
        """
        if df is None or len(df) == 0:
            return df
        out = []
        for chrom, grp in df.groupby("#CHR", sort=False):
            g = grp.sort_values("START")
            starts = g["START"].to_numpy()
            ends = g["END"].to_numpy()
            cs, ce = starts[0], ends[0]
            for s, e in zip(starts[1:], ends[1:]):
                if s <= ce:
                    ce = max(ce, e)
                else:
                    out.append((chrom, cs, ce))
                    cs, ce = s, e
            out.append((chrom, cs, ce))
        return pd.DataFrame(out, columns=["#CHR", "START", "END"])

    del region_df
    if blacklist_df is not None:
        shade_regions(ax, axis, _merge_ranges(blacklist_df), color="gray", alpha=0.2)


def _load_shading(region_bed, blacklist_bed):
    """Read the two optional shading BEDs; either path may be ``None``."""
    return (
        read_bed(region_bed) if region_bed else None,
        read_bed(blacklist_bed) if blacklist_bed else None,
    )


def _suptitle(fig, title, handles=None):
    """Bold page title, placed just above a ``tight_layout``-ed figure.

    *handles* draws a one-row legend on the title's line, at the right edge, so the
    key never covers data.
    """
    fig.subplots_adjust(top=1 - 0.4 / fig.get_figheight())
    y = 1 - 0.12 / fig.get_figheight()
    fig.suptitle(title, fontweight="bold", y=y)
    if handles:
        fig.legend(
            handles=handles,
            loc="center right",
            bbox_to_anchor=(1.0, y),
            ncol=len(handles),
            fontsize=9,
            frameon=False,
            markerscale=1,
        )


def _finish_page(
    fig, title, feature_label, out_file=None, dpi=72, pdf=None, legend_handles=None
):
    """Label, lay out and title one genome-wide page, then write and close it.

    Args:
        fig: The figure to finish.
        title: Bold page super-title.
        feature_label: Feature named in the x-label (``bin``, ``SNP``, ``bb``).
        out_file: Destination when *pdf* is None.
        dpi: Raster resolution.
        pdf: Open ``PdfPages`` to append to; the caller closes it.
        legend_handles: Legend entries for the title line; None draws no legend.
    """
    fig.supxlabel(f"Genome positions (MB) - {feature_label}")
    fig.tight_layout()
    _suptitle(fig, title, handles=legend_handles)
    # NB: bbox_inches keeps artists drawn outside the axes, such as an offset legend
    if pdf is not None:
        pdf.savefig(fig, dpi=dpi, bbox_inches="tight")
    else:
        fig.savefig(out_file, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def _bold_chrnames(ax):
    """Bold the off-axis chromosome-name texts drawn under an ``mb_ticks`` axis."""
    for t in ax.texts:
        t.set_fontweight("bold")


def _hist_with_stats(
    ax, vals, xlabel, header="", ylabel="# segments", clip_q=0.99, sci_x=False
):
    """Histogram of *vals* with a multi-line, mean/median-annotated title.

    *header* is an optional first title line (used for the page-1 panel names); page-2
    rows leave it empty and carry the sample label as a vertical row label instead.
    *sci_x* draws the x-axis in scientific notation (matplotlib's offset multiplier)
    for large-count axes, rather than scaling the values into the label.
    """
    prefix = f"{header}\n" if header else ""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        ax.set_title(f"{prefix}{xlabel}\n(n=0)", fontsize=8, fontweight="bold")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return
    mean, median = float(np.mean(vals)), float(np.median(vals))
    plot_vals = vals
    if clip_q is not None and len(vals) > 1:
        hi = np.quantile(vals, clip_q)
        if hi > 0:
            plot_vals = vals[vals <= hi]
    ax.hist(plot_vals, bins=50, alpha=0.7)
    ax.axvline(mean, color="red", linestyle=":", linewidth=1)
    ax.axvline(median, color="orange", linestyle=":", linewidth=1)
    ax.set_title(
        f"{prefix}{xlabel}\nmean={mean:.1f}, median={median:.1f}",
        fontsize=8,
        fontweight="bold",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if sci_x:
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
