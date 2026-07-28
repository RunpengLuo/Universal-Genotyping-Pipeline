"""Read-depth bias-correction QC plots (rd_correct / count_reads step)."""

import logging

import numpy as np

from scipy.stats import gaussian_kde, pearsonr, spearmanr

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cnplot import adaptive_dot_size, plot_scatter_1d, read_bed

from plot_utils import _get_axis, _shade, _val_full


def _plot_cov_panel(
    ax,
    covariate,
    reads,
    xlabel,
    show_ylabel,
    label,
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
        xs, ys = x[idx], y[idx]
    else:
        xs, ys = x, y

    xlo = xlim[0] if xlim is not None else x.min()
    xhi = xlim[1] if xlim is not None else x.max()
    try:
        kde = gaussian_kde(np.vstack([xs, ys]))
        xgrid = np.linspace(xlo, xhi, 200)
        ygrid = np.linspace(0, ylim * 1.5, 200)
        xx, yy = np.meshgrid(xgrid, ygrid)
        zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
        ax.pcolormesh(xx, yy, zz, shading="gouraud", cmap="Blues", rasterized=True)
        ax.contour(xx, yy, zz, levels=6, colors="steelblue", linewidths=0.5, alpha=0.5)
    except np.linalg.LinAlgError:
        # near-constant covariate -> singular KDE covariance; plain scatter instead
        logging.warning(
            f"KDE failed ({xlabel}, {label}): near-constant covariate; scatter fallback"
        )
        ax.scatter(xs, ys, s=2, color="steelblue", alpha=0.3, rasterized=True)

    mad = np.median(np.abs(y - np.median(y)))
    r_pearson, _ = pearsonr(x, y)
    r_spearman, _ = spearmanr(x, y)
    metrics = f"MAD={mad:.2f}  r={r_pearson:.4f}  rho={r_spearman:.4f}"
    if rmse is not None and is_before:
        metrics = f"RMSE={rmse:.2f}  " + metrics
    logging.info(f"  {xlabel:<8s} {label:<16s}: {metrics}")
    ax.set_xlabel(xlabel)
    if show_ylabel:
        ax.set_ylabel("Observed Readcov")
    ax.set_title(f"{label}\n{metrics}", fontsize=8, fontweight="bold")
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
    labels,
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
    nsamples = len(labels)
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
            for si, label in enumerate(labels):
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
                    label=label,
                    rmse=rmse,
                    is_before=is_before,
                    **kw,
                )
        fig.suptitle(f"{title_prefix}{title}", fontsize=14, fontweight="bold")
        plt.tight_layout()
        pdf.savefig(fig, dpi=150)
        plt.close(fig)
        is_before = False


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
    axis = _get_axis(genome_size, pos_df["#CHR"])
    region_df = read_bed(region_bed) if region_bed else None
    blacklist_df = read_bed(blacklist_bed) if blacklist_bed else None
    s_plot = adaptive_dot_size(len(pos_df), s_base=s)
    alphas = np.full(len(pos_df), alpha)

    for si, label in enumerate(labels):
        fig, axes = plt.subplots(2, 1, figsize=(20, 6), sharex=True)
        for ax, mat, ylim, stage in zip(
            axes,
            [dp_before, dp_after],
            [ylim_before, ylim_after],
            ["before correction", "after correction"],
        ):
            y = mat[:, si] if mat.ndim == 2 else mat
            _shade(ax, axis, region_df, blacklist_df)
            plot_scatter_1d(
                ax,
                pos_df.assign(_y=y),
                axis,
                "_y",
                alphas=alphas,
                markersize=s_plot,
                ylim=(0.0, ylim) if ylim is not None else None,
                ylabel=val_type,
                title=f"{stage} — {_val_full(val_type)}",
                mb_ticks=True,
                show_gaps=False,
            )
            ax.grid(axis="y", alpha=0.2)
        fig.suptitle(str(label), fontsize=12, y=1.0, fontweight="bold")
        fig.supxlabel(f"Genome positions (MB) - {unit}")
        fig.tight_layout()
        pdf.savefig(fig, dpi=dpi)
        plt.close(fig)
