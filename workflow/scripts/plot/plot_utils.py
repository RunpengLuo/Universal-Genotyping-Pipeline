"""Shared plotting helpers: genome axis, region shading, value names, stats histogram.

The common base imported by the per-step plot modules (``plot_genome``,
``plot_count_reads``, ``plot_alleles``, ``plot_combine_counts``,
``plot_genotype_snps``); it creates no figures itself.
"""

import numpy as np
import pandas as pd

from scipy.sparse import issparse

from cnplot import GenomeAxis, shade_regions

from io_utils import get_chr_sizes


def _extract_col(mat, col_idx):
    """Extract a single matrix column as a 1-D numpy array (dense or sparse)."""
    if issparse(mat):
        return np.asarray(mat[:, col_idx].toarray()).ravel()
    return np.asarray(mat[:, col_idx]).ravel()


# Full names for value-type abbreviations, used in figure titles only (ylabels keep
# the short form). AF is the reference-allele frequency (ref / total).
_VAL_TYPE_FULL = {
    "AF": "REF-allele frequency",
    "BAF": "B-allele frequency",
    "RD": "Read depth",
    "RDR": "Read-depth ratio",
}


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
        excluded = [c for c in get_chr_sizes(genome_size) if c not in keep]
        axis = GenomeAxis(
            None, genome_size, excluded_chroms=excluded, collapse_gaps=False
        )
        _AXIS_CACHE[key] = axis
    return axis


def _merge_intervals(df):
    """Union overlapping/touching [START, END) intervals per ``#CHR``.

    Shading each raw interval separately stacks the alpha where they overlap, so
    overlapping masked regions render darker than a single interval. Merging first
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


def _shade(ax, axis, region_df, blacklist_df):
    """Shade masked (blacklist) regions gray; callable regions keep the white background.

    *region_df* (the callable whitelist) is unshaded so the plotted background stays
    white; it is kept in the signature for caller stability. Blacklist intervals are
    unioned first so overlaps do not compound the alpha.
    """
    del region_df
    if blacklist_df is not None:
        shade_regions(ax, axis, _merge_intervals(blacklist_df), color="gray", alpha=0.2)


def _bold_chrnames(ax):
    """Bold the off-axis chromosome-name texts drawn under an ``mb_ticks`` axis."""
    for t in ax.texts:
        t.set_fontweight("bold")


_ASSAY_PLOT_RANK = {"bulkWGS": 0, "bulkWGS-lr": 1, "bulkWES": 2}


def sample_row_order(assays, sample_types, dataset_ids):
    """Row order for stacked sample plots.

    Sorts by assay (WGS < WGS-lr < WES; other assays last), then normal before
    tumor, then ``dataset_id``. Returns the permutation of row indices.
    """

    def key(i):
        return (
            _ASSAY_PLOT_RANK.get(assays[i], len(_ASSAY_PLOT_RANK)),
            0 if sample_types[i] == "normal" else 1,
            str(dataset_ids[i]),
        )

    return sorted(range(len(dataset_ids)), key=key)


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
