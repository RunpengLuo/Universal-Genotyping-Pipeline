"""SNP allele-frequency and depth QC plots (phase_and_concat + binning steps)."""

import os
import logging

import numpy as np
import pandas as pd

from scipy.sparse import issparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

from matrix_utils import dense_col, pseudobulk_by_groups

from plot_genome import plot_1d_sample, plot_1d_multi_sample
from plot_utils import _suptitle


def _af(num, den):
    """``num / den`` as float32, ``NaN`` wherever the denominator is not positive."""
    out = np.full(np.shape(den), np.nan, dtype=np.float32)
    return np.divide(num, den, out=out, where=(den > 0))


def compute_af_per_sample(tot_mtx, b_mtx, i: int):
    """Per-SNP allele frequency for one sample column; ``NaN`` where depth is zero."""
    return _af(dense_col(b_mtx, i), dense_col(tot_mtx, i))


def compute_af_by_groups(tot_mtx, b_mtx, group_idx, n_groups):
    """Per-group pseudobulk allele frequency matrix of shape (n_features, n_groups)."""
    return _af(
        pseudobulk_by_groups(b_mtx, group_idx, n_groups),
        pseudobulk_by_groups(tot_mtx, group_idx, n_groups),
    )


def compute_af_pseudobulk(tot_mtx, b_mtx):
    """Per-SNP allele frequency summed over all cells; ``NaN`` where depth is zero."""

    def row_sum(mat):
        return np.asarray(mat.sum(axis=1)).ravel() if issparse(mat) else mat.sum(axis=1)

    return _af(row_sum(b_mtx), row_sum(tot_mtx))


def plot_snp_depth(
    tot_mtx,
    dataset_ids,
    qc_dir,
    run_id,
    ref_mtx=None,
    b_mtx=None,
    is_bulk=True,
    cell_rep_idx=None,
    name_prefix="",
    pdf: PdfPages | None = None,
    row_order=None,
    sample_id=None,
    max_points=50000,
):
    """Stacked violin+box QC of per-dataset SNP metrics (one axis per metric).

    One row per metric, datasets on a shared x-axis (labels bold, on the last axis
    only): total allele depth (log-y, covered SNPs), then ref-AF and B-AF (linear
    0-1, 0.5 reference line) when *ref_mtx* / *b_mtx* are given. Each dataset is a
    violin with an inner box (median/IQR); values are subsampled to *max_points*
    for the density. The metric name is the (bold) y-label; the page carries a bold
    *sample_id* super-title and no per-axis titles.

    When *pdf* is given the figure is appended to it (and *qc_dir*/*run_id*/
    *name_prefix* are unused); otherwise a standalone PDF is written to *qc_dir*.
    The caller pre-slices the matrices to the SNP rows of interest.

    Parameters
    ----------
    tot_mtx, ref_mtx, b_mtx : ndarray or sparse
        Total depth, ref-allele, and B-allele count matrices (SNPs x samples/cells).
        *ref_mtx* / *b_mtx* are optional; each adds its own violin row.
    dataset_ids : list[str]
        Sample / replicate identifiers (violin x-order).
    is_bulk : bool
        If True, matrix columns are samples. If False, columns are cells; see
        *cell_rep_idx*.
    cell_rep_idx : np.ndarray or None
        Length-n_cells rep index per cell, consulted only when ``is_bulk=False``:
        cells are pseudobulked within each rep, else all cells collapse to one.
    row_order : list[int] or None
        Dataset permutation applied before drawing (match the 1-D scatter order).
    sample_id : str or None
        Sample/patient id for the bold page super-title.
    max_points : int
        Per-dataset subsample cap for the violin density.
    """
    logging.info("QC analysis - plot SNP depth violin")

    def _resolve(mat):
        """Per-dataset matrix in dataset (column) space, or None."""
        if mat is None:
            return None
        if is_bulk:
            return mat
        if cell_rep_idx is not None:
            return pseudobulk_by_groups(mat, cell_rep_idx, len(dataset_ids))
        return (
            np.asarray(mat.sum(axis=1))
            if issparse(mat)
            else mat.sum(axis=1, keepdims=True)
        )

    depth_mat = _resolve(tot_mtx)
    ref_count_mat = _resolve(ref_mtx)
    b_count_mat = _resolve(b_mtx)
    if is_bulk or cell_rep_idx is not None:
        labels = list(dataset_ids)
    else:
        labels = ["pseudobulk"]

    if row_order is not None and len(row_order) == len(labels):
        labels = [labels[i] for i in row_order]
        depth_mat = depth_mat[:, row_order]
        if ref_count_mat is not None:
            ref_count_mat = ref_count_mat[:, row_order]
        if b_count_mat is not None:
            b_count_mat = b_count_mat[:, row_order]

    rng = np.random.default_rng(0)

    def _sub(v):
        return rng.choice(v, max_points, replace=False) if len(v) > max_points else v

    def _long(value_fn):
        rows = []
        for ci, label in enumerate(labels):
            rows.append(pd.DataFrame({"dataset": label, "val": _sub(value_fn(ci))}))
        return pd.concat(rows, ignore_index=True)

    def _depth(ci):
        total = dense_col(depth_mat, ci).astype(np.float64)
        return total[total > 0]

    def _frac(count_mat):
        def fn(ci):
            total = dense_col(depth_mat, ci).astype(np.float64)
            cov = total > 0
            return dense_col(count_mat, ci).astype(np.float64)[cov] / total[cov]

        return fn

    # (ylabel, kind, long_df); kind "log" => log depth axis, "frac" => 0-1 with 0.5 line
    metrics = [("Total allele depth", "log", _long(_depth))]
    if ref_count_mat is not None:
        metrics.append(("Ref allele frequency", "frac", _long(_frac(ref_count_mat))))
    if b_count_mat is not None:
        metrics.append(("B-allele frequency", "frac", _long(_frac(b_count_mat))))

    n = len(metrics)
    fig, axes = plt.subplots(
        n, 1, figsize=(max(6.0, 1.3 * len(labels)), 3.5 * n), sharex=True, squeeze=False
    )
    axes = axes[:, 0]
    for ax, (ylabel, kind, df) in zip(axes, metrics):
        sns.violinplot(
            data=df, x="dataset", y="val", order=labels, inner="box", cut=0, ax=ax
        )
        ax.set_xlabel("")
        ax.set_ylabel(ylabel, fontweight="bold")
        if kind == "log":
            ax.set_yscale("log")
        else:
            ax.axhline(0.5, color="grey", linestyle="--", linewidth=0.8)
            ax.set_ylim(0, 1)
    plt.setp(axes[-1].get_xticklabels(), rotation=45, ha="right", fontweight="bold")

    fig.tight_layout()
    if sample_id:
        _suptitle(fig, sample_id)
    if pdf is not None:
        pdf.savefig(fig)
        plt.close(fig)
        logging.info("added SNP depth violin page")
        return
    stem = f"{name_prefix}.snp_depth" if name_prefix else "snp_depth"
    out_path = os.path.join(qc_dir, f"{stem}.{run_id}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    logging.info(f"saved SNP depth violin to {out_path}")


def plot_allele_freqs(
    pos_df,
    dataset_ids,
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
    sample_id=None,
    row_order=None,
):
    """Generate genome-wide allele-frequency scatter plots.

    Output mode is chosen by the combination of ``apply_pseudobulk`` and
    ``cell_rep_idx``:

    - ``apply_pseudobulk=False`` — columns of the matrices are samples;
      multi-row scatter, one row per ``dataset_ids`` entry.
    - ``apply_pseudobulk=True`` and ``cell_rep_idx`` provided — cells are
      pseudobulked within each rep; multi-row scatter, one row per rep.
    - ``apply_pseudobulk=True`` and ``cell_rep_idx`` is ``None`` — all cells
      collapse into a single pseudobulk page.

    Parameters
    ----------
    pos_df : pd.DataFrame
        SNP/bin position DataFrame with ``#CHR`` and ``POS`` (or ``START``/``END``).
    dataset_ids : list[str]
        Replicate identifiers; used as row labels.
    tot_mtx, b_mtx : sparse or ndarray
        Total depth and B-allele count matrices.
    cell_rep_idx : np.ndarray or None
        Length-n_cells int array mapping each cell column to a rep index in
        ``dataset_ids``. See behavior matrix above.
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
            sample_id=sample_id,
            region_bed=region_bed,
            blacklist_bed=blacklist_bed,
            pdf=pdf,
        )
        return

    if per_rep_pseudobulk:
        af_mat = compute_af_by_groups(tot_mtx, b_mtx, cell_rep_idx, len(dataset_ids))
    else:
        _tot_mtx = tot_mtx.tocsc() if issparse(tot_mtx) else tot_mtx
        _b_mtx = b_mtx.tocsc() if issparse(b_mtx) else b_mtx
        af_mat = np.column_stack(
            [
                compute_af_per_sample(_tot_mtx, _b_mtx, i)
                for i in range(len(dataset_ids))
            ]
        )
    stem = f"af_{allele}_{unit}{suffix}"
    stem = f"{name_prefix}.{stem}" if name_prefix else stem
    plot_file = os.path.join(plot_dir, f"{stem}.{run_id}.pdf")
    plot_1d_multi_sample(
        pos_df,
        af_mat,
        list(dataset_ids),
        genome_size,
        plot_file,
        unit=unit,
        val_type=val_type,
        sample_id=sample_id,
        row_order=row_order,
        region_bed=region_bed,
        blacklist_bed=blacklist_bed,
        pdf=pdf,
    )
    return
