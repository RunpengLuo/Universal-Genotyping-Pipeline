"""Genome-wide 1D scatter primitives (single- and multi-sample), on the cnplot axis."""

import logging

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

from cnplot import adaptive_dot_size, plot_scatter_1d, read_bed

from plot_utils import _bold_chrnames, _get_axis, _shade, _val_full


def plot_1d_multi_sample(
    pos_df: pd.DataFrame,
    mat: np.ndarray,
    labels: list,
    genome_size: str,
    out_file: str,
    unit="window",
    val_type="RD",
    s=4,
    dpi=72,
    alpha=0.6,
    min_ylim=0.0,
    max_ylim=None,
    sample_id: str | None = None,
    row_order: list | None = None,
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
    row_order : list[int] or None
        Row permutation applied to *labels* and *mat* columns before drawing.
    """
    if row_order is not None:
        labels = [labels[i] for i in row_order]
        mat = mat[:, row_order]
    n_samples = len(labels)
    logging.info(
        f"genome-wide {unit}-level {val_type} multi-sample plot "
        f"({n_samples} samples), out_file={out_file}"
    )
    axis = _get_axis(genome_size, pos_df["#CHR"])
    region_df = read_bed(region_bed) if region_bed else None
    blacklist_df = read_bed(blacklist_bed) if blacklist_bed else None
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
            ylabel=label,
            plot_chrname=is_last,
            mb_ticks=True,
            show_gaps=False,
        )
        if is_last:
            _bold_chrnames(ax)

    val_name = _val_full(val_type)
    title = f"{sample_id} - {val_name}" if sample_id else val_name
    fig.supxlabel(f"Genome positions (MB) - {unit}")
    fig.tight_layout()
    fig.subplots_adjust(top=1 - 0.4 / fig.get_figheight())
    fig.suptitle(title, fontweight="bold", y=1 - 0.12 / fig.get_figheight())
    if pdf is not None:
        pdf.savefig(fig, dpi=dpi)
    else:
        fig.savefig(out_file, dpi=dpi)
    plt.close(fig)


def plot_1d_sample(
    pos_df: pd.DataFrame,
    val: np.ndarray,
    genome_size: str,
    out_file: str,
    unit="SNP",
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
    """
    logging.info(f"genome-wide {unit}-level {val_type} plot, out_file={out_file}")
    axis = _get_axis(genome_size, pos_df["#CHR"])
    region_df = read_bed(region_bed) if region_bed else None
    blacklist_df = read_bed(blacklist_bed) if blacklist_bed else None

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

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        hue = np.where(mask, mask_labels[0], mask_labels[1])
        palette = {mask_labels[0]: mask_colors[0], mask_labels[1]: mask_colors[1]}
        plot_scatter_1d(
            ax,
            pos_df.assign(_y=val, _hue=hue),
            axis,
            "_y",
            hue="_hue",
            palette=palette,
            alphas=np.full(len(pos_df), 0.8),
            markersize=s_plot,
            href=0.5 if is_frac else None,
            ylim=ylim,
            ylabel=val_type,
            mb_ticks=True,
            show_gaps=False,
        )
        on, off = int((m & mask).sum()), int((m & ~mask).sum())
        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color=mask_colors[0],
                label=f"{mask_labels[0]} ({on})",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="",
                color=mask_colors[1],
                label=f"{mask_labels[1]} ({off})",
            ),
        ]
        ax.legend(handles=handles, loc="upper right", fontsize=8, markerscale=1)
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
    fig.supxlabel(f"Genome positions (MB) - {unit}")
    fig.tight_layout()
    fig.subplots_adjust(top=1 - 0.4 / fig.get_figheight())
    fig.suptitle(title, fontweight="bold", y=1 - 0.12 / fig.get_figheight())
    if pdf is not None:
        pdf.savefig(fig, dpi=dpi)
    else:
        fig.savefig(out_file, dpi=dpi)
    plt.close(fig)
    return
