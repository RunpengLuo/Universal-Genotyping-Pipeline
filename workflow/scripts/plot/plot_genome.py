"""Genome-wide 1D scatter primitives, shared by the per-step plot modules.

Last update: 2026-08-08

Functions:
- plot_1d_sample: one genome-wide track, optionally split by a mask
- plot_1d_multi_sample: one row per sample on a shared axis
"""

import logging

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

from cnplot import adaptive_dot_size, plot_scatter_1d

from plot_utils import (
    _bold_chrnames,
    _finish_page,
    _get_axis,
    _load_shading,
    _shade,
    _val_full,
)


def plot_1d_multi_sample(
    pos_df: pd.DataFrame,
    mat: np.ndarray,
    labels: list,
    genome_size: str,
    out_file: str,
    feature_label="bin",
    val_type="RD",
    s=4,
    dpi=72,
    alpha=0.6,
    min_ylim=0.0,
    max_ylim=None,
    sample_id: str | None = None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
):
    """Multi-sample 1-D genome-wide scatter plot: all chromosomes on one page,
    one row per sample.

    The value name is a bold page super-title (``{sample_id} - {value name}`` when
    *sample_id* is given). Chromosome names are drawn (bold) under the last row only.

    Parameters
    ----------
    mat : np.ndarray
        (n_windows, n_samples) value matrix.
    labels : list[str]
        Sample labels, length == mat.shape[1].
    sample_id : str or None
        Sample/patient id for the page super-title.
    """
    n_samples = len(labels)
    logging.info(
        f"genome-wide {feature_label}-level {val_type} multi-sample plot "
        f"({n_samples} samples), out_file={out_file}"
    )
    axis = _get_axis(genome_size, pos_df["#CHR"])
    region_df, blacklist_df = _load_shading(region_bed, blacklist_bed)
    s_plot = adaptive_dot_size(len(pos_df), s_base=s)
    alphas = np.full(len(pos_df), alpha)
    is_frac = val_type in ("AF", "BAF")

    fig, axes = plt.subplots(
        nrows=n_samples,
        ncols=1,
        figsize=(20, 3 * n_samples),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]

    for si, (ax, label) in enumerate(zip(axes, labels)):
        is_last = si == n_samples - 1
        y = mat[:, si] if mat.ndim == 2 else mat
        _shade(ax, axis, region_df, blacklist_df)
        plot_scatter_1d(
            ax,
            pos_df.assign(_y=y),
            axis,
            "_y",
            alphas=alphas,
            markersize=s_plot,
            href=0.5 if is_frac else None,
            ylim=(-0.05, 1.05)
            if is_frac
            else ((min_ylim, max_ylim) if max_ylim is not None else None),
            # rotated axis: break "{dataset_id} {assay} {T|N}" after the dataset_id
            ylabel=label.replace(" ", "\n", 1),
            plot_chrname=is_last,
            mb_ticks=True,
            show_gaps=False,
        )
        if is_last:
            _bold_chrnames(ax)

    val_name = _val_full(val_type)
    title = f"{sample_id} - {val_name}" if sample_id else val_name
    _finish_page(fig, title, feature_label, out_file=out_file, dpi=dpi, pdf=pdf)


def plot_1d_sample(
    pos_df: pd.DataFrame,
    val: np.ndarray,
    genome_size: str,
    out_file: str,
    feature_label="SNP",
    val_type="BAF",
    s=4,
    dpi=72,
    alpha=0.6,
    figsize=(20, 3),
    min_ylim=0.0,
    max_ylim=None,
    mask: np.ndarray | None = None,
    mask_labels=("kept", "filtered"),
    mask_colors=("blue", "red"),
    groups: np.ndarray | None = None,
    group_colors: dict | None = None,
    sample_id: str | None = None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
):
    """Single-sample 1-D genome-wide scatter plot: all chromosomes on one page.

    The value name is a bold page super-title (``{sample_id} - {value name}`` when
    *sample_id* is given). Chromosome names are drawn bold under the axis.

    When *mask* is given, ``mask``-true points use ``mask_colors[0]``/``mask_labels[0]``
    and ``mask``-false points use ``mask_colors[1]``/``mask_labels[1]``, with a legend.
    *groups* generalizes that to any number of classes: a per-point label array plus a
    ``{label: color}`` map, whose insertion order is the legend order. Either way the
    legend carries each class's point count, so classes may share a color and stay
    countable.
    """
    logging.info(
        f"genome-wide {feature_label}-level {val_type} plot, out_file={out_file}"
    )
    axis = _get_axis(genome_size, pos_df["#CHR"])
    region_df, blacklist_df = _load_shading(region_bed, blacklist_bed)

    m = np.isfinite(val)
    s_plot = adaptive_dot_size(int(m.sum()), s_base=s)
    is_frac = val_type in ("AF", "BAF")
    ylim = (
        (-0.05, 1.05)
        if is_frac
        else ((min_ylim, max_ylim) if max_ylim is not None else None)
    )

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    _shade(ax, axis, region_df, blacklist_df)

    if groups is None and mask is not None:
        mask = np.asarray(mask, dtype=bool)
        groups = np.where(mask, mask_labels[0], mask_labels[1])
        group_colors = {mask_labels[0]: mask_colors[0], mask_labels[1]: mask_colors[1]}

    if groups is not None:
        groups = np.asarray(groups)
        plot_scatter_1d(
            ax,
            pos_df.assign(_y=val, _hue=groups),
            axis,
            "_y",
            hue="_hue",
            palette=dict(group_colors),
            alphas=np.full(len(pos_df), 0.8),
            markersize=s_plot,
            href=0.5 if is_frac else None,
            ylim=ylim,
            ylabel=val_type,
            mb_ticks=True,
            show_gaps=False,
        )
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color=color,
                label=f"{label} ({int((m & (groups == label)).sum())})",
            )
            for label, color in group_colors.items()
        ]
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(1.005, 1.0),
            borderaxespad=0,
            fontsize=8,
            markerscale=1,
            frameon=False,
        )
    else:
        plot_scatter_1d(
            ax,
            pos_df.assign(_y=val),
            axis,
            "_y",
            alphas=np.full(len(pos_df), alpha),
            markersize=s_plot,
            href=0.5 if is_frac else None,
            ylim=ylim,
            ylabel=val_type,
            mb_ticks=True,
            show_gaps=False,
        )
    _bold_chrnames(ax)

    val_name = _val_full(val_type)
    title = f"{sample_id} - {val_name}" if sample_id else val_name
    _finish_page(fig, title, feature_label, out_file=out_file, dpi=dpi, pdf=pdf)
    return
