"""Bin-level combine_counts QC: segmentation, genome-wide RDR/BAF, RDR-vs-BAF.

Last update: 2026-08-28

Functions:
- plot_loh_density: het-SNP density per tile and the clonal-LOH regions called from it
- plot_segmentation_qc: bb length and per-dataset raw count histograms
- plot_rdr_baf: one page per tumor, depth then RDR then BAF
- plot_rdr_baf_2d: RDR-vs-BAF cloud with marginal densities
"""

import logging

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.lines import Line2D

from cnplot import adaptive_dot_size, plot_scatter_1d, plot_scatter_2d

from segmentation_utils import dense_observation

from plot_utils import (
    _bold_chrnames,
    _finish_page,
    _get_axis,
    _hist_with_stats,
    _load_shading,
    _observation_labels,
    _shade,
)


def plot_segmentation_qc(
    seg_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    x_count_mat,
    rd_count_mat,
    b_count_mat,
    tot_count_mat,
    out_file: str | None = None,
    pdf: PdfPages | None = None,
    sample_id: str = "",
    dpi: int = 150,
    is_loh=None,
):
    """Two-page segmentation QC histograms for combine_counts output.

    Page 1 — segment length (kbp) over all segments, split non-LOH / clonal-LOH when
      *is_loh* is given, since the two are binned on different criteria.
    Page 2 — one row per dataset_id, four histograms of raw counts: native counts,
      read starts, B-allele counts, total-allele counts. The count axes use scientific
      notation (matplotlib's offset multiplier) rather than a scaled axis label.
      Each row is labelled ``{dataset_id}\\n{assay_type} {T|N}`` on the rotated row axis;
      the patient id is the page super-title.

    Parameters
    ----------
    seg_df : pd.DataFrame
        Segmentation table with ``#CHR``, ``START``, ``END`` (bb.tsv.gz schema).
    sample_df : pd.DataFrame
        One row per count-matrix column (per dataset_id), with ``dataset_id``, ``assay_type``
        and ``sample_type``. Row order must match the columns of the count matrices.
    x_count_mat, rd_count_mat, b_count_mat, tot_count_mat : ndarray or sparse, (n_seg, n_datasets)
        Native, read-start, B-allele, and total-allele counts per segment per dataset_id;
        columns aligned to *sample_df* rows.
    out_file, pdf : see the other ``plot_*`` functions. Exactly one is used.
    is_loh : (n_seg,) bool array or None
        Clonal-LOH flag per segment; None draws one length panel instead of two.
    """
    logging.info("QC analysis - plot segmentation QC histograms")

    n_datasets = len(sample_df)

    _own_pdf = pdf is None
    pdf_pages = PdfPages(out_file) if _own_pdf else pdf

    # ---- page 1: segment length ----
    lengths_kbp = (seg_df["END"].to_numpy() - seg_df["START"].to_numpy()) / 1000.0

    if is_loh is None:
        panels = [(lengths_kbp, "Segment length")]
    else:
        is_loh = np.asarray(is_loh, dtype=bool)
        panels = [
            (
                lengths_kbp[~is_loh],
                f"Segment length, non-LOH (n={int((~is_loh).sum())})",
            ),
            (
                lengths_kbp[is_loh],
                f"Segment length, clonal LOH (n={int(is_loh.sum())})",
            ),
        ]
    fig1, axes1 = plt.subplots(
        1, len(panels), figsize=(5.5 * len(panels), 4), squeeze=False
    )
    for ax1, (vals, header) in zip(axes1[0], panels):
        _hist_with_stats(ax1, vals, "segment length (kbp)", header)
    fig1.suptitle(
        f"Segmentation QC — {len(seg_df)} segments", fontsize=11, fontweight="bold"
    )
    fig1.tight_layout()
    pdf_pages.savefig(fig1, dpi=dpi)
    plt.close(fig1)

    # ---- page 2: per-dataset_id count histograms ----
    fig2, axes = plt.subplots(
        nrows=max(n_datasets, 1),
        ncols=4,
        figsize=(20, 3 * max(n_datasets, 1)),
        squeeze=False,
    )
    for ri in range(n_datasets):
        row = sample_df.iloc[ri]
        # shown once per row as a bold vertical "row super-title"
        row_label = (
            f"{row.get('dataset_id', '')}\n{row.get('assay_type', '')} "
            f"{str(row.get('sample_type', ''))[:1].upper()}"
        )
        _hist_with_stats(
            axes[ri, 0], dense_observation(x_count_mat, ri), "aligned bases", sci_x=True
        )
        _hist_with_stats(
            axes[ri, 1],
            dense_observation(rd_count_mat, ri),
            "read-start count",
            sci_x=True,
        )
        _hist_with_stats(
            axes[ri, 2],
            dense_observation(b_count_mat, ri),
            "B-allele count",
            sci_x=True,
        )
        _hist_with_stats(
            axes[ri, 3],
            dense_observation(tot_count_mat, ri),
            "total allele count",
            sci_x=True,
        )
        axes[ri, 0].annotate(
            row_label,
            xy=(0, 0.5),
            xytext=(-axes[ri, 0].yaxis.labelpad - 22, 0),
            xycoords=axes[ri, 0].yaxis.label,
            textcoords="offset points",
            ha="right",
            va="center",
            rotation=90,
            fontweight="bold",
            fontsize=9,
        )
    fig2.suptitle(
        f"{sample_id} - per-dataset counts" if sample_id else "per-dataset counts",
        fontsize=11,
        fontweight="bold",
    )
    fig2.tight_layout()
    fig2.subplots_adjust(left=0.18)
    pdf_pages.savefig(fig2, dpi=dpi)
    plt.close(fig2)

    if _own_pdf:
        pdf_pages.close()
        logging.info(f"saved segmentation QC histograms to {out_file}")


# padded so a bb stamped at BAF 0 sits inside the axes, not on the spine; the ticks stay
# on the [0, 1] the quantity is defined over
BAF_XLIM = (-0.01, 1.01)
BAF_XTICKS = np.linspace(0, 1, 6)


def _rdr_ylim(rdr_mat):
    """Upper RDR axis limit: the 99th percentile rounded to an integer, +1, +10%."""
    return (np.round(np.nanquantile(rdr_mat, 0.99)).astype(int) + 1) * 1.1


def plot_rdr_baf(
    pos_df: pd.DataFrame,
    rdr_mat: np.ndarray,
    baf_mat: np.ndarray,
    depth_tumor_mat: np.ndarray,
    depth_normal_mat: np.ndarray,
    sample_id: str,
    dataset_ids: list,
    assay_types: list,
    rdr_base_dataset_ids: list,
    genome_size: str,
    out_file: str,
    feature_label="bb",
    s=4,
    dpi=72,
    alpha=0.7,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
    tumor_color="tab:red",
    normal_color="tab:blue",
    is_loh=None,
    loh_color="tab:green",
):
    """Genome-wide read-depth + RDR + BAF plot: one page per tumor, three rows.

    Each page stacks, top to bottom: matched-normal & tumor read depth on a shared
    axis (normal drawn first, tumor overlaid), tumor RDR, tumor BAF. Chromosome
    names are drawn (bold) under the BAF row only; the page carries a bold
    super-title instead of per-axis titles, with the point key on the same line.

    Args:
        pos_df: Position DataFrame with ``#CHR`` and ``START``/``END`` columns.
        rdr_mat: (N, T) tumor RDR per bin.
        baf_mat: (N, T) tumor BAF per bin.
        depth_tumor_mat: (N, T) tumor read depth per bin.
        depth_normal_mat: (N, T) matched-normal read depth per bin; an all-NaN
            column means the tumor has no matched normal and only tumor depth is
            drawn.
        sample_id: Sample/patient id, opening every page super-title.
        dataset_ids: Tumor dataset id per column, length T.
        assay_types: Tumor assay type per column, length T.
        rdr_base_dataset_ids: The dataset id each tumor's RDR is divided by, or
            ``None`` for a median-normalized tumor; length T. Sets the ``/ base`` half
            of the page title.
        genome_size: Path to chromosome sizes file.
        out_file: Output PDF path; used only when ``pdf`` is ``None``.
        feature_label: Feature named in the x-label (e.g. ``"bb"``).
        region_bed: Path to whitelist BED for background shading.
        blacklist_bed: Path to blacklist BED for background shading.
        pdf: External ``PdfPages``; pages are appended and the caller closes it.
        tumor_color: Tumor dot color.
        normal_color: Normal dot color.
        is_loh: (N,) bool array marking clonal-LOH bins, drawn in *loh_color*. Their
            BAF is stamped at 0 rather than measured, so it must read as a separate
            class; None draws every bin as a tumor bin.
        loh_color: Clonal-LOH dot color.
    """
    n_tumors = len(dataset_ids)
    logging.info(
        f"genome-wide {feature_label}-level depth+RDR+BAF plot "
        f"({n_tumors} tumors), out_file={out_file}"
    )
    labels = _observation_labels(dataset_ids, assay_types, ["tumor"] * n_tumors)
    rdr_ylim = _rdr_ylim(rdr_mat)
    axis = _get_axis(genome_size, pos_df["#CHR"])
    region_df, blacklist_df = _load_shading(region_bed, blacklist_bed)
    s_plot = adaptive_dot_size(len(pos_df), s_base=s)
    alphas = np.full(len(pos_df), alpha)
    if is_loh is None:
        point_colors = tumor_color
        n_loh = 0
    else:
        is_loh = np.asarray(is_loh, dtype=bool)
        point_colors = np.where(is_loh, loh_color, tumor_color)
        n_loh = int(is_loh.sum())

    _own_pdf = pdf is None
    pdf_pages = PdfPages(out_file) if _own_pdf else pdf
    for si in range(n_tumors):
        base = rdr_base_dataset_ids[si]
        title = f"{sample_id} - {labels[si]}"
        if base:
            title = f"{title} / {base}"
        fig, (ax_dp, ax_rdr, ax_baf) = plt.subplots(
            nrows=3, ncols=1, figsize=(20, 9), sharex=True
        )

        # --- read-depth row: matched normal first, tumor overlaid ---
        dp_normal = depth_normal_mat[:, si]
        has_normal = np.isfinite(dp_normal).any()
        _shade(ax_dp, axis, region_df, blacklist_df)
        if has_normal:
            plot_scatter_1d(
                ax_dp,
                pos_df.assign(_y=dp_normal),
                axis,
                "_y",
                colors=normal_color,
                alphas=alphas,
                markersize=s_plot,
                ylabel="Read-depth",
                plot_chrname=False,
                mb_ticks=True,
                show_gaps=False,
            )
        plot_scatter_1d(
            ax_dp,
            pos_df.assign(_y=depth_tumor_mat[:, si]),
            axis,
            "_y",
            colors=point_colors,
            alphas=alphas,
            markersize=s_plot,
            ylabel="Read-depth",
            plot_chrname=False,
            mb_ticks=True,
            show_gaps=False,
        )
        # the page key, drawn on the title line by _finish_page rather than over the data
        handles = [
            Line2D([0], [0], marker="o", linestyle="", color=tumor_color, label="tumor")
        ]
        if has_normal:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    color=normal_color,
                    label="normal",
                )
            )
        if n_loh:
            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    color=loh_color,
                    label="clonal LOH",
                )
            )

        # --- RDR row: tumor only ---
        _shade(ax_rdr, axis, region_df, blacklist_df)
        plot_scatter_1d(
            ax_rdr,
            pos_df.assign(_y=rdr_mat[:, si]),
            axis,
            "_y",
            colors=point_colors,
            alphas=alphas,
            markersize=s_plot,
            ylim=(0, rdr_ylim),
            ylabel="RDR",
            plot_chrname=False,
            mb_ticks=True,
            show_gaps=False,
        )

        # --- BAF row: tumor only, chr names bold on this last axis ---
        _shade(ax_baf, axis, region_df, blacklist_df)
        plot_scatter_1d(
            ax_baf,
            pos_df.assign(_y=baf_mat[:, si]),
            axis,
            "_y",
            colors=point_colors,
            alphas=alphas,
            markersize=s_plot,
            href=0.5,
            ylim=(-0.05, 1.05),
            ylabel="BAF",
            mb_ticks=True,
            show_gaps=False,
        )
        _bold_chrnames(ax_baf)

        _finish_page(
            fig, title, feature_label, dpi=dpi, pdf=pdf_pages, legend_handles=handles
        )
    if _own_pdf:
        pdf_pages.close()


def plot_rdr_baf_2d(
    rdr_mat,
    baf_mat,
    sample_id: str,
    dataset_ids: list,
    assay_types: list,
    out_file: str | None = None,
    dpi: int = 150,
    pdf: PdfPages | None = None,
):
    """RDR-vs-BAF joint scatter, one page per tumor (no copy-number landmarks).

    A diagnostic companion to the genome-wide ``plot_rdr_baf``: the observed
    per-bin cloud with marginal densities, via ``cnplot.plot_scatter_2d``. The BAF axis
    is padded to :data:`BAF_XLIM` while its ticks stay on ``[0, 1]``, so a bin at 0 or 1
    is drawn inside the axes rather than half under the spine.

    Args:
        rdr_mat, baf_mat: (n_bins, T) RDR and BAF per bin per tumor.
        sample_id: Sample/patient id, opening every page title.
        dataset_ids: Tumor dataset id per column, length T.
        assay_types: Tumor assay type per column, length T.
        out_file: Output PDF path; used only when ``pdf`` is None.
        dpi: Raster resolution.
        pdf: External PdfPages; pages are appended and the caller closes it.
    """
    labels = _observation_labels(dataset_ids, assay_types, ["tumor"] * len(dataset_ids))
    logging.info(f"QC analysis - RDR-vs-BAF 2D scatter ({len(labels)} tumors)")
    rdr_ylim = _rdr_ylim(rdr_mat)
    _own_pdf = pdf is None
    pdf_pages = PdfPages(out_file) if _own_pdf else pdf
    for si, label in enumerate(labels):
        rdr = rdr_mat[:, si] if rdr_mat.ndim == 2 else rdr_mat
        baf = baf_mat[:, si] if baf_mat.ndim == 2 else baf_mat
        finite = np.isfinite(rdr) & np.isfinite(baf)
        obs = pd.DataFrame({"RDR": rdr[finite], "BAF": baf[finite]})
        grid = plot_scatter_2d(
            obs,
            "BAF",
            "RDR",
            xlim=BAF_XLIM,
            ylim=(0, rdr_ylim),
            refline_x=0.5,
            refline_y=1.0,
            xlabel="BAF",
            ylabel="RDR",
            title=f"{sample_id} - {label}",
        )
        grid.ax_joint.set_xticks(BAF_XTICKS)
        pdf_pages.savefig(grid.figure, dpi=dpi)
        plt.close(grid.figure)
    if _own_pdf:
        pdf_pages.close()


def plot_loh_density(
    tiles: pd.DataFrame,
    rate_ref: float,
    rate_loh: float,
    loh_df: pd.DataFrame,
    genome_size: str,
    sample_id: str,
    out_file: str | None = None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
    dpi: int = 72,
    neutral_color="tab:red",
    loh_color="tab:green",
):
    """Genome-wide het-SNP density per tile, coloured by the state the chain gave it.

    One page. The two fitted rates are drawn as horizontal references, so a reader can
    see whether the decode followed the data or the prior: a real clonal-LOH region sits
    on the lower line with nothing in between.

    Args:
        tiles: Tile frame from ``loh_utils.build_loh_tiles`` plus an ``is_loh`` column.
        rate_ref, rate_loh: Fitted neutral and assumed LOH het rate, per Mb.
        loh_df: The reported intervals, for the count in the title.
        genome_size: Path to chromosome sizes file.
        sample_id: Sample id, opening the page title.
        out_file: Output PDF path; used only when *pdf* is None.
        region_bed, blacklist_bed: Background shading.
        pdf: External ``PdfPages``; pages are appended and the caller closes it.
        dpi: Raster resolution.
        neutral_color, loh_color: Dot colour per state.
    """
    exposure_mb = tiles["exposure_bp"].to_numpy() / 1e6
    with np.errstate(invalid="ignore", divide="ignore"):
        density = np.where(
            exposure_mb > 0, tiles["n_het"].to_numpy() / exposure_mb, np.nan
        )
    is_loh = tiles["is_loh"].to_numpy(dtype=bool)
    logging.info(
        f"QC analysis - clonal-LOH density over {len(tiles)} tiles, "
        f"{int(is_loh.sum())} called"
    )

    axis = _get_axis(genome_size, tiles["#CHR"])
    region_df, blacklist_df = _load_shading(region_bed, blacklist_bed)
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 4))
    _shade(ax, axis, region_df, blacklist_df)
    plot_scatter_1d(
        ax,
        tiles.assign(_y=density),
        axis,
        "_y",
        colors=np.where(is_loh, loh_color, neutral_color),
        alphas=np.full(len(tiles), 0.7),
        markersize=adaptive_dot_size(len(tiles), s_base=6),
        ylim=(-0.05 * rate_ref, 1.6 * rate_ref),
        ylabel="het SNPs / Mb",
        mb_ticks=True,
        show_gaps=False,
    )
    for rate, style, label in (
        (rate_ref, "--", f"neutral {rate_ref:.0f}/Mb"),
        (rate_loh, ":", f"LOH {rate_loh:.1f}/Mb"),
    ):
        ax.axhline(rate, linestyle=style, linewidth=1, color="black", zorder=1)
        ax.annotate(
            label,
            xy=(1.0, rate),
            xycoords=("axes fraction", "data"),
            fontsize=8,
            ha="right",
            va="bottom",
        )
    _bold_chrnames(ax)
    handles = [
        Line2D([0], [0], marker="o", ls="", color=neutral_color, label="neutral tile"),
        Line2D([0], [0], marker="o", ls="", color=loh_color, label="clonal-LOH tile"),
    ]
    span_mb = (
        float((loh_df["END"] - loh_df["START"]).sum()) / 1e6 if len(loh_df) else 0.0
    )
    _finish_page(
        fig,
        f"{sample_id} - clonal LOH: {len(loh_df)} regions, {span_mb:.0f} Mbp",
        "tile",
        out_file=out_file,
        dpi=dpi,
        pdf=pdf,
        legend_handles=handles,
    )
