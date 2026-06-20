"""Plotting utilities for QC and diagnostic visualizations."""

import os
import logging

import numpy as np
import pandas as pd

from scipy.stats import gaussian_kde, pearsonr, spearmanr
from scipy.sparse import issparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from io_utils import get_chr_sizes, read_region_file
from utils import adaptive_dot_size
from combine_counts_utils import (
    compute_af_pseudobulk,
    compute_af_per_sample,
    compute_af_by_groups,
    pseudobulk_by_groups,
)


# ---------------------------------------------------------------------------
# rd_correct_utils plots
# ---------------------------------------------------------------------------


def _plot_cov_panel(
    ax,
    covariate,
    reads,
    xlabel,
    show_ylabel,
    rep_id,
    rmse=None,
    is_before=False,
    xlim=None,
    xticks=None,
):
    """KDE density scatter of readcov vs a covariate on a single axes."""
    valid = (reads > 0) & np.isfinite(covariate)
    if valid.sum() < 20:
        ax.set_visible(False)
        return
    x, y = covariate[valid], reads[valid]
    ylim = np.nanquantile(y, 0.99) * 1.1

    n_pts = len(x)
    rng = np.random.default_rng(0)
    if n_pts > 20000:
        idx = rng.choice(n_pts, size=20000, replace=False)
        kde = gaussian_kde(np.vstack([x[idx], y[idx]]))
    else:
        kde = gaussian_kde(np.vstack([x, y]))

    xlo = xlim[0] if xlim is not None else x.min()
    xhi = xlim[1] if xlim is not None else x.max()
    xgrid = np.linspace(xlo, xhi, 200)
    ygrid = np.linspace(0, ylim * 1.5, 200)
    xx, yy = np.meshgrid(xgrid, ygrid)
    positions = np.vstack([xx.ravel(), yy.ravel()])
    zz = kde(positions).reshape(xx.shape)

    ax.pcolormesh(xx, yy, zz, shading="gouraud", cmap="Blues", rasterized=True)
    ax.contour(xx, yy, zz, levels=6, colors="steelblue", linewidths=0.5, alpha=0.5)

    mad = np.median(np.abs(y - np.median(y)))
    r_pearson, _ = pearsonr(x, y)
    r_spearman, _ = spearmanr(x, y)
    metrics = f"MAD={mad:.2f}  r={r_pearson:.4f}  rho={r_spearman:.4f}"
    if rmse is not None and is_before:
        metrics = f"RMSE={rmse:.2f}  " + metrics
    logging.info(f"  {xlabel:<8s} {rep_id:<8s}: {metrics}")
    ax.set_xlabel(xlabel)
    if show_ylabel:
        ax.set_ylabel("Observed Readcov")
    ax.set_title(f"{rep_id}\n{metrics}", fontsize=8)
    if xlim is not None:
        ax.set_xlim(*xlim)
    else:
        ax.set_xlim(x.min(), x.max())
    if xticks is not None:
        ax.set_xticks(xticks)
    ax.set_ylim(0, ylim * 1.1)


def plot_rd_2d_kde(
    gc,
    dp_before,
    dp_after,
    rep_ids,
    pdf,
    gc_rmse=None,
    mappability=None,
    repliseq=None,
    title_prefix="",
):
    """Two-page PDF: before/after correction KDE density plots.

    Each page has up to 3 rows (GC, MAP, RT) x nsamples columns.

    Parameters
    ----------
    gc_rmse : list of float or None
        Per-sample RMSE from the GC fit. Shown on the "Before" panel only.
    mappability : np.ndarray or None
        Per-window mappability values. If provided, a MAP row is added.
    repliseq : np.ndarray or None
        Per-window replication timing values. If provided, an RT row is added.
    title_prefix : str
        Optional prefix for page titles (e.g. ``"target — "``).
    """
    nsamples = len(rep_ids)
    panel_w = max(5, 5 * nsamples)

    # Build list of (row_label, covariate, xlabel)
    rows = [("GC", gc, "GC Content")]
    if mappability is not None:
        rows.append(("MAP", mappability, "Mappability"))
    if repliseq is not None:
        rows.append(("RT", repliseq, "Replication Timing"))
    nrows = len(rows)

    is_before = True
    for title, dp_mat in [
        ("Before Correction", dp_before),
        ("After Correction", dp_after),
    ]:
        logging.info(f"====={title}=====")
        fig, axes = plt.subplots(
            nrows,
            nsamples,
            figsize=(panel_w, 5 * nrows),
            squeeze=False,
        )
        for ri, (row_label, covariate, xlabel) in enumerate(rows):
            for si, rep_id in enumerate(rep_ids):
                rmse = (
                    gc_rmse[si] if (gc_rmse is not None and row_label == "GC") else None
                )
                kw = {}
                if row_label == "MAP":
                    kw = {"xlim": (-0.2, 1.2), "xticks": np.arange(0, 1.1, 0.2)}
                _plot_cov_panel(
                    axes[ri, si],
                    covariate,
                    dp_mat[:, si],
                    xlabel,
                    show_ylabel=(si == 0),
                    rep_id=rep_id,
                    rmse=rmse,
                    is_before=is_before,
                    **kw,
                )
        fig.suptitle(f"{title_prefix}{title}", fontsize=14)
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        is_before = False


# ---------------------------------------------------------------------------
# count_reads_utils plots
# ---------------------------------------------------------------------------


def plot_rd_1d_scatter(
    pos_df,
    dp_before,
    dp_after,
    labels,
    genome_size,
    pdf,
    unit="window",
    val_type="RD",
    ylim_before=None,
    ylim_after=None,
    s=4,
    dpi=72,
    alpha=0.6,
    region_bed=None,
    blacklist_bed=None,
):
    """One page per sample: top = before correction, bottom = after correction."""
    genome_x, chrom_bounds, total_len = _genome_coords(pos_df, genome_size)
    region_by_chr = _parse_bed_by_chr(region_bed)
    blacklist_by_chr = _parse_bed_by_chr(blacklist_bed)
    s_plot = adaptive_dot_size(len(pos_df), s_base=s)

    for si, label in enumerate(labels):
        fig, axes = plt.subplots(2, 1, figsize=(20, 6), sharex=True)
        for ai, (ax, mat, ylim, title) in enumerate(zip(
            axes,
            [dp_before, dp_after],
            [ylim_before, ylim_after],
            ["before correction", "after correction"],
        )):
            y = mat[:, si] if mat.ndim == 2 else mat
            m = np.isfinite(y)
            _add_chrom_decorations(ax, chrom_bounds, total_len, region_by_chr, blacklist_by_chr)
            if m.any():
                ax.scatter(genome_x[m], y[m], s=s_plot, alpha=alpha, rasterized=True)
            if ylim is not None:
                ax.set_ylim(0.0, ylim)
            ax.grid(axis="y", alpha=0.2)
            ax.set_title(f"{title} — {val_type} ({unit})", fontsize=10)
            ax.set_ylabel(val_type, fontsize=10)
        fig.suptitle(f"SAMPLE {label}", fontsize=12, y=1.0)
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)


def _genome_coords(pos_df, genome_size):
    """Compute genome-wide x-coordinates, chromosome offsets, and boundaries.

    Returns (genome_x, chrom_bounds, total_len) where chrom_bounds is a list
    of (chrom_label, offset, chrom_size) tuples.
    """
    chrom_sizes = get_chr_sizes(genome_size)
    ch = pos_df["#CHR"].to_numpy()
    if "POS" in pos_df.columns:
        pos = pos_df["POS"].to_numpy()
    else:
        pos = ((pos_df["START"].to_numpy() + pos_df["END"].to_numpy()) // 2).astype(
            np.int64
        )

    # Determine chromosome order from data
    seen = []
    for c in ch:
        if len(seen) == 0 or seen[-1] != c:
            seen.append(c)

    offsets = {}
    cum = 0
    chrom_bounds = []
    for c in seen:
        sz = chrom_sizes.get(c)
        if sz is None:
            continue
        offsets[c] = cum
        chrom_bounds.append((c, cum, sz))
        cum += sz

    genome_x = np.array(
        [offsets.get(c, 0) + p for c, p in zip(ch, pos)], dtype=np.float64
    )
    return genome_x, chrom_bounds, cum


def _add_chrom_decorations(
    ax, chrom_bounds, total_len, region_by_chr, blacklist_by_chr
):
    """Add chromosome boundaries, labels, and region/blacklist shading to an axis."""
    for chrom, offset, sz in chrom_bounds:
        ax.axvline(offset, color="black", linewidth=0.5, alpha=0.3)
        if chrom in region_by_chr:
            for s, e in region_by_chr[chrom]:
                ax.axvspan(
                    offset + s, offset + e, color="lightblue", alpha=0.15, zorder=0
                )
        if chrom in blacklist_by_chr:
            for s, e in blacklist_by_chr[chrom]:
                ax.axvspan(
                    offset + s, offset + e, color="lightcoral", alpha=0.2, zorder=0
                )
        label = chrom.replace("chr", "")
        ax.text(
            offset + sz / 2,
            -0.06,
            label,
            ha="center",
            va="top",
            fontsize=7,
            transform=ax.get_xaxis_transform(),
        )
    ax.set_xlim(0, total_len)
    ax.set_xticks([])


def _parse_bed_by_chr(bed_path):
    """Read a BED file and return {chrom: [(start, end), ...]}."""
    result = {}
    if bed_path is not None:
        df = read_region_file(bed_path)
        for ch, grp in df.groupby("#CHR", sort=False):
            result[ch] = list(zip(grp["START"], grp["END"]))
    return result


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
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
):
    """Multi-sample 1-D genome-wide scatter plot: all chromosomes on one page,
    one row per sample.

    Parameters
    ----------
    mat : np.ndarray
        (n_windows, n_samples) value matrix.
    labels : list[str]
        Sample labels, length == mat.shape[1].
    """
    n_samples = len(labels)
    logging.info(
        f"genome-wide {unit}-level {val_type} multi-sample plot "
        f"({n_samples} samples), out_file={out_file}"
    )
    genome_x, chrom_bounds, total_len = _genome_coords(pos_df, genome_size)
    region_by_chr = _parse_bed_by_chr(region_bed)
    blacklist_by_chr = _parse_bed_by_chr(blacklist_bed)

    s_plot = adaptive_dot_size(len(pos_df), s_base=s)

    fig, axes = plt.subplots(
        nrows=n_samples,
        ncols=1,
        figsize=(20, 3 * n_samples),
        sharex=True,
        squeeze=False,
    )
    axes = axes[:, 0]

    for si, (ax, label) in enumerate(zip(axes, labels)):
        y = mat[:, si] if mat.ndim == 2 else mat
        m = np.isfinite(y)

        _add_chrom_decorations(
            ax, chrom_bounds, total_len, region_by_chr, blacklist_by_chr
        )

        if m.any():
            ax.scatter(genome_x[m], y[m], s=s_plot, alpha=alpha, rasterized=True)

        if val_type in ["AF", "BAF"]:
            ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
            ax.set_ylim(-0.05, 1.05)
        elif max_ylim is not None:
            ax.set_ylim(min_ylim, max_ylim)

        ax.set_ylabel(label, fontsize=10)
        if val_type not in ["AF", "BAF"]:
            ax.grid(axis="y", alpha=0.2)

        if si == 0:
            ax.set_title(f"{val_type} ({unit})")

    fig.tight_layout()
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
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
):
    """Single-sample 1-D genome-wide scatter plot: all chromosomes on one page."""
    logging.info(f"genome-wide {unit}-level {val_type} plot, out_file={out_file}")
    genome_x, chrom_bounds, total_len = _genome_coords(pos_df, genome_size)
    region_by_chr = _parse_bed_by_chr(region_bed)
    blacklist_by_chr = _parse_bed_by_chr(blacklist_bed)

    m = np.isfinite(val)
    s_plot = adaptive_dot_size(int(m.sum()), s_base=s)

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    _add_chrom_decorations(ax, chrom_bounds, total_len, region_by_chr, blacklist_by_chr)

    if mask is not None:
        kept = m & mask
        filt = m & ~mask
        if filt.any():
            ax.scatter(
                genome_x[filt],
                val[filt],
                s=s_plot,
                alpha=0.8,
                color="red",
                rasterized=True,
                label=f"filtered ({filt.sum()})",
            )
        if kept.any():
            ax.scatter(
                genome_x[kept],
                val[kept],
                s=s_plot,
                alpha=0.8,
                color="blue",
                rasterized=True,
                label=f"kept ({kept.sum()})",
            )
        ax.legend(loc="upper right", fontsize=8, markerscale=2)
    else:
        if m.any():
            ax.scatter(genome_x[m], val[m], s=s_plot, alpha=alpha, rasterized=True)

    if val_type in ["AF", "BAF"]:
        ax.axhline(0.5, color="grey", linestyle=":", linewidth=1)
        ax.set_ylim(-0.05, 1.05)
    elif max_ylim is not None:
        ax.set_ylim(min_ylim, max_ylim)

    ax.set_ylabel(val_type)
    if val_type not in ["AF", "BAF"]:
        ax.grid(axis="y", alpha=0.2)
    ax.set_title(f"{val_type} ({unit})")
    fig.tight_layout()
    if pdf is not None:
        pdf.savefig(fig, dpi=dpi)
    else:
        fig.savefig(out_file, dpi=dpi)
    plt.close(fig)
    return


def plot_rdr_baf(
    pos_df: pd.DataFrame,
    rdr_mat: np.ndarray,
    baf_mat: np.ndarray,
    labels: list,
    genome_size: str,
    out_file: str,
    unit="bb",
    s=4,
    dpi=72,
    alpha=0.6,
    rdr_ylim=None,
    region_bed: str | None = None,
    blacklist_bed: str | None = None,
    pdf: PdfPages | None = None,
):
    """Combined RDR + BAF genome-wide plot: one page per sample, two rows.

    Produces one figure page per sample with RDR on top and BAF on the bottom.

    Args:
        pos_df: Position DataFrame with ``#CHR`` and ``START``/``END`` columns.
        rdr_mat: Array of shape (N, S) with RDR values per bin per sample.
        baf_mat: Array of shape (N, S) with BAF values per bin per sample.
        labels: Sample labels, length S.
        genome_size: Path to chromosome sizes file.
        out_file: Output PDF path. Used only when ``pdf`` is ``None``; ignored
            when an external ``PdfPages`` object is supplied via ``pdf``.
        unit: Feature unit label for axis titles (e.g. ``"bb"``).
        rdr_ylim: Upper y-axis limit for the RDR panel. No limit if ``None``.
        region_bed: Path to whitelist BED for background shading.
        blacklist_bed: Path to blacklist BED for background shading.
        pdf: External ``PdfPages`` object. When provided, pages are appended to
            it and the caller is responsible for closing the object. When
            ``None``, a new PDF is created at ``out_file`` and closed on return.
    """
    n_samples = len(labels)
    logging.info(
        f"genome-wide {unit}-level RDR+BAF plot "
        f"({n_samples} samples), out_file={out_file}"
    )
    genome_x, chrom_bounds, total_len = _genome_coords(pos_df, genome_size)
    region_by_chr = _parse_bed_by_chr(region_bed)
    blacklist_by_chr = _parse_bed_by_chr(blacklist_bed)
    s_plot = adaptive_dot_size(len(pos_df), s_base=s)

    _own_pdf = pdf is None
    pdf_pages = PdfPages(out_file) if _own_pdf else pdf
    for si, label in enumerate(labels):
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(20, 6), sharex=True)

        # RDR
        ax_rdr = axes[0]
        y_rdr = rdr_mat[:, si] if rdr_mat.ndim == 2 else rdr_mat
        m_rdr = np.isfinite(y_rdr)
        _add_chrom_decorations(
            ax_rdr, chrom_bounds, total_len, region_by_chr, blacklist_by_chr
        )
        if m_rdr.any():
            ax_rdr.scatter(
                genome_x[m_rdr], y_rdr[m_rdr], s=s_plot, alpha=alpha, rasterized=True
            )
        if rdr_ylim is not None:
            ax_rdr.set_ylim(0, rdr_ylim)
        ax_rdr.set_ylabel("RDR")
        ax_rdr.grid(axis="y", alpha=0.2)
        ax_rdr.set_title(f"{label} — RDR + BAF ({unit})")

        # BAF
        ax_baf = axes[1]
        y_baf = baf_mat[:, si] if baf_mat.ndim == 2 else baf_mat
        m_baf = np.isfinite(y_baf)
        _add_chrom_decorations(
            ax_baf, chrom_bounds, total_len, region_by_chr, blacklist_by_chr
        )
        if m_baf.any():
            ax_baf.scatter(
                genome_x[m_baf], y_baf[m_baf], s=s_plot, alpha=alpha, rasterized=True
            )
        ax_baf.axhline(0.5, color="grey", linestyle=":", linewidth=1)
        ax_baf.set_ylim(-0.05, 1.05)
        ax_baf.set_ylabel("BAF")

        fig.tight_layout()
        pdf_pages.savefig(fig, dpi=dpi)
        plt.close(fig)
    if _own_pdf:
        pdf_pages.close()


# ---------------------------------------------------------------------------
# combine_counts_utils plots
# ---------------------------------------------------------------------------


# native per-segment count label by assay (ax-1 of the segmentation QC page-2).
# Bulk passes total aligned bases (sum of per-window depth*length), shown in Mbp;
# ASSAY_COUNT_SCALE divides the raw values to match the label's unit.
ASSAY_COUNT_LABEL = {
    "bulkWGS": "DNA aligned bases (Mbp)",
    "bulkWGS-lr": "DNA aligned bases (Mbp)",
    "bulkWES": "DNA aligned bases (Mbp)",
    "scRNA": "UMI count",
    "VISIUM": "UMI count",
    "VISIUM3prime": "UMI count",
    "scATAC": "ATAC fragment count",
}
ASSAY_COUNT_SCALE = {
    "bulkWGS": 1e6,
    "bulkWGS-lr": 1e6,
    "bulkWES": 1e6,
}


def _seg_gene_counts(seg_df, gene_count, gene_col):
    """Resolve a per-segment gene count from an explicit array or a seg_df column."""
    if gene_count is not None:
        return np.asarray(gene_count, dtype=float)
    if "n_genes" in seg_df.columns:
        return seg_df["n_genes"].to_numpy(dtype=float)
    if gene_col in seg_df.columns:
        def _count(v):
            if not isinstance(v, str) or v == "":
                return 0
            return sum(1 for g in v.split(";") if g and g != "intergenic")
        return seg_df[gene_col].map(_count).to_numpy(dtype=float)
    return None


def _hist_with_stats(ax, vals, xlabel, header="", ylabel="# segments", clip_q=0.99):
    """Histogram of *vals* with a multi-line, mean/median-annotated title.

    *header* is an optional first title line (used for the page-1 panel names); page-2
    rows leave it empty and carry the sample label as a vertical row label instead.
    """
    prefix = f"{header}\n" if header else ""
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        ax.set_title(f"{prefix}{xlabel}\n(n=0)", fontsize=8)
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
    ax.set_title(f"{prefix}{xlabel}\nmean={mean:.1f}, median={median:.1f}", fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def plot_segmentation_qc(
    seg_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    x_count_mat,
    b_count_mat,
    tot_count_mat,
    out_file: str | None = None,
    pdf: PdfPages | None = None,
    gene_count=None,
    gene_col: str = "feature_id",
    dpi: int = 150,
):
    """Two-page segmentation QC histograms for combine_counts output.

    Page 1 — two histograms over all segments:
      (i) segment length (kbp), (ii) per-segment gene count.
    Page 2 — one row per REP_ID, three histograms:
      native counts (DNA reads / UMI / ATAC fragments, by assay), B-allele counts,
      total-allele counts. Every subplot title is multi-line and carries
      ``{SAMPLE} {REP_ID} {assay_type} {sample_type}`` plus mean/median.

    Parameters
    ----------
    seg_df : pd.DataFrame
        Segmentation table with ``#CHR``, ``START``, ``END`` (bb.tsv.gz schema).
        For the gene-count panel it may carry a ``;``-joined ``gene_col`` or a numeric
        ``n_genes`` column; otherwise pass *gene_count* explicitly.
    sample_df : pd.DataFrame
        One row per count-matrix column (per REP_ID), with ``REP_ID``, ``assay_type``,
        ``sample_type`` and a sample-name column (``SAMPLE_NAME`` or ``SAMPLE``). Row
        order must match the columns of the count matrices.
    x_count_mat, b_count_mat, tot_count_mat : ndarray or sparse, (n_seg, n_rep)
        Native, B-allele, and total-allele counts per segment per rep; columns aligned
        to *sample_df* rows.
    gene_count : array-like (n_seg,) or None
        Optional explicit per-segment gene count (overrides derivation from *seg_df*).
    out_file, pdf : see the other ``plot_*`` functions. Exactly one is used.
    """
    logging.info("QC analysis - plot segmentation QC histograms")

    def _col(mat, j):
        if issparse(mat):
            return np.asarray(mat[:, j].toarray()).ravel()
        return np.asarray(mat[:, j]).ravel()

    name_col = "SAMPLE_NAME" if "SAMPLE_NAME" in sample_df.columns else "SAMPLE"
    n_rep = len(sample_df)

    _own_pdf = pdf is None
    pdf_pages = PdfPages(out_file) if _own_pdf else pdf

    # ---- page 1: segment length + per-segment gene count ----
    lengths_kbp = (seg_df["END"].to_numpy() - seg_df["START"].to_numpy()) / 1000.0
    genes_per_seg = _seg_gene_counts(seg_df, gene_count, gene_col)
    fig1, ax1 = plt.subplots(1, 2, figsize=(11, 4))
    _hist_with_stats(ax1[0], lengths_kbp, "segment length (kbp)", "Segment length")
    if genes_per_seg is not None:
        _hist_with_stats(ax1[1], genes_per_seg, "# genes / segment", "Genes per segment")
    else:
        ax1[1].set_title("Genes per segment\n(no gene annotation)", fontsize=8)
        ax1[1].set_xlabel("# genes / segment")
        ax1[1].set_ylabel("# segments")
    fig1.suptitle(f"Segmentation QC — {len(seg_df)} segments", fontsize=11)
    fig1.tight_layout()
    pdf_pages.savefig(fig1, dpi=dpi)
    plt.close(fig1)

    # ---- page 2: per-rep count histograms ----
    fig2, axes = plt.subplots(
        nrows=max(n_rep, 1), ncols=3, figsize=(15, 3 * max(n_rep, 1)), squeeze=False
    )
    for ri in range(n_rep):
        row = sample_df.iloc[ri]
        assay = str(row.get("assay_type", ""))
        # 2-line sample label, shown once per row as a bold vertical "row super-title"
        row_label = (
            f"{row.get(name_col, '')} {row.get('REP_ID', '')}\n"
            f"{assay} {row.get('sample_type', '')}"
        )
        x_label = ASSAY_COUNT_LABEL.get(assay, "count")
        x_scale = ASSAY_COUNT_SCALE.get(assay, 1.0)
        _hist_with_stats(axes[ri, 0], _col(x_count_mat, ri) / x_scale, x_label)
        _hist_with_stats(axes[ri, 1], _col(b_count_mat, ri), "B-allele count")
        _hist_with_stats(axes[ri, 2], _col(tot_count_mat, ri), "total allele count")
        axes[ri, 0].annotate(
            row_label,
            xy=(0, 0.5), xytext=(-axes[ri, 0].yaxis.labelpad - 22, 0),
            xycoords=axes[ri, 0].yaxis.label, textcoords="offset points",
            ha="right", va="center", rotation=90, fontweight="bold", fontsize=9,
        )
    fig2.tight_layout()
    fig2.subplots_adjust(left=0.18)
    pdf_pages.savefig(fig2, dpi=dpi)
    plt.close(fig2)

    if _own_pdf:
        pdf_pages.close()
        logging.info(f"saved segmentation QC histograms to {out_file}")


def plot_snp_depth_histogram(
    tot_mtx,
    rep_ids,
    qc_dir,
    run_id,
    ref_mtx=None,
    is_bulk=True,
    cell_rep_idx=None,
    name_prefix="",
):
    """Plot per-sample histograms of total allele depth and ref-AF at SNP positions.

    The caller is responsible for pre-slicing *tot_mtx* (and *ref_mtx*) to the
    desired SNP rows before calling this function.

    Parameters
    ----------
    tot_mtx : ndarray or sparse
        Total depth matrix (SNPs x samples/cells), already filtered to the
        SNPs of interest.
    rep_ids : list[str]
        Sample / replicate identifiers.
    qc_dir : str
        Output directory for the PDF.
    run_id : str
        Run identifier appended to the output filename before the extension.
    ref_mtx : ndarray or sparse, optional
        Ref-allele count matrix (same shape as *tot_mtx*). When provided, a
        second column of ref-AF histograms is added to the figure.
    is_bulk : bool
        If True, columns of the matrices are samples — one row per sample.
        If False, columns are cells; see *cell_rep_idx*.
    cell_rep_idx : np.ndarray or None
        Length-n_cells int array mapping each cell column to a rep index in
        ``rep_ids``. Only consulted when ``is_bulk=False``. When provided,
        cells are pseudobulked within each rep — one row per rep. When
        ``None``, all cells collapse into a single "pseudobulk" row.
    """
    logging.info("QC analysis - plot SNP depth histogram")

    has_af = ref_mtx is not None
    ncols = 2 if has_af else 1

    if is_bulk:
        depth_mat = tot_mtx
        ref_count_mat = ref_mtx
        row_labels = list(rep_ids)
    elif cell_rep_idx is not None:
        n_groups = len(rep_ids)
        depth_mat = pseudobulk_by_groups(tot_mtx, cell_rep_idx, n_groups)
        ref_count_mat = (
            pseudobulk_by_groups(ref_mtx, cell_rep_idx, n_groups)
            if has_af
            else None
        )
        row_labels = list(rep_ids)
    else:
        # legacy global-pseudobulk path: one row aggregating all cells
        if issparse(tot_mtx):
            depth_mat = np.asarray(tot_mtx.sum(axis=1))
        else:
            depth_mat = tot_mtx.sum(axis=1, keepdims=True)
        if has_af:
            if issparse(ref_mtx):
                ref_count_mat = np.asarray(ref_mtx.sum(axis=1))
            else:
                ref_count_mat = ref_mtx.sum(axis=1, keepdims=True)
        else:
            ref_count_mat = None
        row_labels = ["pseudobulk"]
    nrows = len(row_labels)

    fig, axes = plt.subplots(
        nrows=nrows, ncols=ncols, figsize=(5 * ncols, 3 * nrows), squeeze=False
    )

    def _stats_title(name, depth_vals):
        if len(depth_vals) == 0:
            return f"{name} (n=0)"
        return (
            f"{name}\n"
            f"mean={depth_vals.mean():.1f}, median={np.median(depth_vals):.1f}, "
            f"min={depth_vals.min():.0f}, max={depth_vals.max():.0f}"
        )

    def _get_col(mat, col_idx):
        """Extract a single column as a 1-D numpy array."""
        if issparse(mat):
            return np.asarray(mat[:, col_idx].toarray()).ravel()
        return np.asarray(mat[:, col_idx]).ravel()

    for ri, label in enumerate(row_labels):
        # --- depth histogram ---
        ax_depth = axes[ri, 0]
        depth = _get_col(depth_mat, ri)
        if len(depth) > 0:
            clip_threshold = np.percentile(depth, 99)
            ax_depth.hist(depth[depth <= clip_threshold], bins=50, alpha=0.7)
        ax_depth.set_title(_stats_title(label, depth), fontsize=9)
        ax_depth.set_xlabel("Total allele depth")
        ax_depth.set_ylabel("# SNPs")

        # --- ref-AF histogram ---
        if has_af:
            ax_af = axes[ri, 1]
            total_depth = _get_col(depth_mat, ri).astype(np.float64)
            ref_depth = _get_col(ref_count_mat, ri).astype(np.float64)
            covered_mask = total_depth > 0
            ref_af = np.full(len(total_depth), np.nan)
            ref_af[covered_mask] = ref_depth[covered_mask] / total_depth[covered_mask]
            ref_af_valid = ref_af[~np.isnan(ref_af)]
            if len(ref_af_valid) > 0:
                ax_af.hist(ref_af_valid, bins=50, range=(0, 1), alpha=0.7)
            ax_af.axvline(0.5, color="red", linestyle="--", linewidth=0.8)
            ax_af.set_xlim(0, 1)
            ax_af.set_title(f"{label} ref-AF", fontsize=9)
            ax_af.set_xlabel("Ref allele frequency")
            ax_af.set_ylabel("# SNPs")

    fig.tight_layout()
    stem = f"{name_prefix}.snp_depth_hist" if name_prefix else "snp_depth_hist"
    out_path = os.path.join(qc_dir, f"{stem}.{run_id}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    logging.info(f"saved SNP depth histogram to {out_path}")


def plot_allele_freqs(
    pos_df,
    rep_ids,
    tot_mtx,
    b_mtx,
    genome_size,
    plot_dir,
    apply_pseudobulk=False,
    allele="ref",
    unit="SNP",
    suffix="",
    snp_mask=None,
    region_bed=None,
    blacklist_bed=None,
    run_id="",
    pdf: PdfPages | None = None,
    cell_rep_idx=None,
    name_prefix="",
):
    """Generate genome-wide allele-frequency scatter plots.

    Output mode is chosen by the combination of ``apply_pseudobulk`` and
    ``cell_rep_idx``:

    - ``apply_pseudobulk=False`` — columns of the matrices are samples;
      multi-row scatter, one row per ``rep_ids`` entry.
    - ``apply_pseudobulk=True`` and ``cell_rep_idx`` provided — cells are
      pseudobulked within each rep; multi-row scatter, one row per rep.
    - ``apply_pseudobulk=True`` and ``cell_rep_idx`` is ``None`` — all cells
      collapse into a single pseudobulk page.

    Parameters
    ----------
    pos_df : pd.DataFrame
        SNP/bin position DataFrame with ``#CHR`` and ``POS`` (or ``START``/``END``).
    rep_ids : list[str]
        Replicate identifiers; used as row labels.
    tot_mtx, b_mtx : sparse or ndarray
        Total depth and B-allele count matrices.
    cell_rep_idx : np.ndarray or None
        Length-n_cells int array mapping each cell column to a rep index in
        ``rep_ids``. See behavior matrix above.
    genome_size : str
        Path to chromosome sizes file.
    plot_dir : str
        Output directory for PDF plots.
    allele : str
        Allele label for filenames (e.g., ``"ref"``, ``"B"``).
    unit : str
        Feature unit label (e.g., ``"SNP"``, ``"bin"``).
    suffix : str
        Optional filename suffix.
    region_bed, blacklist_bed : str or None
        Optional BED files for background shading.
    """
    per_rep_pseudobulk = apply_pseudobulk and cell_rep_idx is not None
    # B allele is phased -> label as BAF; ref allele is unphased -> AF.
    val_type = "BAF" if allele == "B" else "AF"
    logging.info(
        f"QC analysis - plot {allele}-{unit} allele frequency, "
        f"apply_pseudobulk={apply_pseudobulk}, per_rep_pseudobulk={per_rep_pseudobulk}"
    )

    if apply_pseudobulk and not per_rep_pseudobulk:
        af = compute_af_pseudobulk(tot_mtx, b_mtx)
        stem = f"af_{allele}_{unit}.pseudobulk{suffix}"
        stem = f"{name_prefix}.{stem}" if name_prefix else stem
        plot_file = os.path.join(plot_dir, f"{stem}.{run_id}.pdf")
        plot_1d_sample(
            pos_df,
            af,
            genome_size,
            plot_file,
            unit=unit,
            val_type=val_type,
            mask=snp_mask,
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
            pdf=pdf,
        )
        return

    if per_rep_pseudobulk:
        af_mat = compute_af_by_groups(tot_mtx, b_mtx, cell_rep_idx, len(rep_ids))
    else:
        _tot_mtx = tot_mtx.tocsc() if issparse(tot_mtx) else tot_mtx
        _b_mtx = b_mtx.tocsc() if issparse(b_mtx) else b_mtx
        af_mat = np.column_stack(
            [compute_af_per_sample(_tot_mtx, _b_mtx, i) for i in range(len(rep_ids))]
        )
    stem = f"af_{allele}_{unit}{suffix}"
    stem = f"{name_prefix}.{stem}" if name_prefix else stem
    plot_file = os.path.join(plot_dir, f"{stem}.{run_id}.pdf")
    plot_1d_multi_sample(
        pos_df,
        af_mat,
        list(rep_ids),
        genome_size,
        plot_file,
        unit=unit,
        val_type=val_type,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        pdf=pdf,
    )
    return
