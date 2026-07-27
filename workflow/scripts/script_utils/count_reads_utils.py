"""Utility functions for window/bin depth processing.

Used by rd_correct.py and combine_counts.py.
"""

import numpy as np


def compute_gc_rd_stats(mat, gc_vals, labels, n_gc_bins=100):
    """Compute per-label Pearson/Spearman corr(RD, GC) and std of binned median RD.

    Parameters
    ----------
    mat : np.ndarray
        (n_bins, n_samples) depth matrix.
    gc_vals : np.ndarray
        Per-bin GC fraction (same length as mat rows).
    labels : list[str]
        Column labels (sample/rep IDs).
    n_gc_bins : int
        Number of equal-width GC bins in [0, 1].

    Returns
    -------
    gc_corr : dict[str, tuple[float, float]]
        {label: (pearson_r, spearman_r)}
    gc_bin_median_std : dict[str, float]
        {label: std of per-GC-bin median RD (A_GC)}
    """
    from scipy.stats import pearsonr, spearmanr

    gc_bins = np.linspace(0, 1, n_gc_bins + 1)
    gc_corr = {}
    gc_bin_median_std = {}

    for i, label in enumerate(labels):
        v = mat[:, i] if mat.ndim == 2 else mat
        valid = np.isfinite(v) & np.isfinite(gc_vals)

        if valid.sum() > 2:
            pr_val, _ = pearsonr(gc_vals[valid], v[valid])
            sr_val, _ = spearmanr(gc_vals[valid], v[valid])
            gc_corr[label] = (pr_val, sr_val)
        else:
            gc_corr[label] = (np.nan, np.nan)

        bin_idx = np.digitize(gc_vals, gc_bins) - 1
        bin_idx = np.clip(bin_idx, 0, n_gc_bins - 1)
        medians = []
        for b in range(n_gc_bins):
            mask = (bin_idx == b) & np.isfinite(v)
            if mask.any():
                medians.append(np.median(v[mask]))
        gc_bin_median_std[label] = float(np.std(medians)) if medians else np.nan

    return gc_corr, gc_bin_median_std
