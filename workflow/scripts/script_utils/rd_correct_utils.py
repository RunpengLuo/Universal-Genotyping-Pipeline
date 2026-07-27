"""Utility functions for per-window read depth bias correction.

Contains HMMcopy-style LOWESS correction.
"""

import logging

import numpy as np
from scipy.interpolate import interp1d
from statsmodels.nonparametric.smoothers_lowess import lowess


def correct_readcount_lowess(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    samplesize=50000,
    routlier=0.01,
    doutlier=0.001,
    min_mappability=0.9,
    lowess_frac_tight=0.03,
    lowess_frac_smooth=0.3,
    grid_size=1001,
    seed=42,
):
    """HMMcopy-style correctReadcount: GC + optional mappability + optional
    replication timing stages (ACEseq-style sequential LOWESS correction).

    Parameters
    ----------
    reads : np.ndarray
        Raw read counts per window (1-D).
    gc : np.ndarray
        Per-window GC fraction in [0, 1].
    mappability : np.ndarray or None
        Per-window mappability values in [0, 1]. If None, stage 2 is skipped.
    repliseq : np.ndarray or None
        Per-window consensus replication timing score. If None, stage 3 is
        skipped.  Higher values = earlier replication = higher expected
        coverage in cycling cells.
    samplesize : int
        Max number of ideal bins for loess fitting.
    routlier : float
        Upper quantile for read-count outlier removal.
    doutlier : float
        Quantile for GC/mappability domain outlier removal (top/bottom).
    min_mappability : float
        Minimum mappability for ideal bins.
    lowess_frac_tight : float
        Bandwidth fraction for the first (tight) LOWESS pass.
    lowess_frac_smooth : float
        Bandwidth fraction for the second (smoothing) LOWESS pass.
    grid_size : int
        Number of grid points for LOWESS interpolation.
    seed : int
        Base random seed for subsampling ideal bins (incremented per stage).

    Returns
    -------
    np.ndarray
        Corrected read counts (same length as *reads*).  NaN where
        correction is not possible.
    float
        RMSE from the GC LOWESS fit (root mean squared error between raw
        depth and LOWESS prediction at valid bins).
    """
    reads = reads.astype(np.float64)
    n = len(reads)

    def _dedup_sorted(xy):
        """Remove duplicate x-values from sorted LOWESS output to avoid
        division-by-zero in interp1d."""
        _, idx = np.unique(xy[:, 0], return_index=True)
        return xy[idx]

    def _dedup_input(y, x):
        """Average y-values for duplicate x-values to avoid degenerate
        local regressions (division by zero) inside LOWESS."""
        ux, inv = np.unique(x, return_inverse=True)
        if len(ux) == len(x):
            return y, x
        uy = np.zeros(len(ux), dtype=np.float64)
        np.add.at(uy, inv, y)
        counts = np.bincount(inv).astype(np.float64)
        uy /= counts
        return uy, ux

    def _fit_lowess_interp(y, x, grid):
        """Tight LOWESS -> grid smooth -> final interpolator."""
        y, x = _dedup_input(y, x)
        if len(x) < 2:
            return None
        s1 = lowess(y, x, frac=lowess_frac_tight, return_sorted=True)
        s1 = _dedup_sorted(s1)
        if len(s1) < 2:
            return None
        s1_fn = interp1d(
            s1[:, 0],
            s1[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        s2 = lowess(s1_fn(grid), grid, frac=lowess_frac_smooth, return_sorted=True)
        s2 = _dedup_sorted(s2)
        if len(s2) < 2:
            return None
        return interp1d(
            s2[:, 0],
            s2[:, 1],
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )

    def _apply_stage(
        prev,
        covariate,
        cov_lo,
        cov_hi,
        grid,
        seed,
        stage_name="",
        extra_ideal_mask=None,
    ):
        """Apply one LOWESS correction stage and return (corrected, rmse)."""
        valid = (prev > 0) & np.isfinite(covariate)
        val_hi = np.nanquantile(prev[valid], 1.0 - routlier)
        ideal = valid & (prev <= val_hi) & (covariate >= cov_lo) & (covariate <= cov_hi)
        if extra_ideal_mask is not None:
            ideal &= extra_ideal_mask

        ideal_idx = np.where(ideal)[0]
        if ideal_idx.size < 10:
            logging.warning(
                f"    {stage_name}: only {ideal_idx.size} ideal bins; skipping correction"
            )
            return prev.copy(), 0.0
        if ideal_idx.size > samplesize:
            rng = np.random.default_rng(seed)
            ideal_idx = rng.choice(ideal_idx, size=samplesize, replace=False)

        interp_fn = _fit_lowess_interp(prev[ideal_idx], covariate[ideal_idx], grid)
        if interp_fn is None:
            logging.warning(
                f"    {stage_name}: LOWESS fit failed (too few distinct values); skipping correction"
            )
            return prev.copy(), 0.0
        predicted = interp_fn(covariate)

        # RMSE between raw values and LOWESS prediction (valid bins only)
        valid_pred = (predicted > 0) & (prev > 0)
        n_valid_pred = int(valid_pred.sum())
        rmse = (
            float(np.sqrt(np.mean((prev[valid_pred] - predicted[valid_pred]) ** 2)))
            if n_valid_pred > 0
            else 0.0
        )

        with np.errstate(invalid="ignore", divide="ignore"):
            corrected = np.where(
                (predicted > 0) & (prev > 0),
                prev / predicted,
                np.nan,
            )
        valid_corr = (predicted > 0) & (prev > 0)
        if valid_corr.any():
            scale = np.median(prev[valid_corr]) / np.median(corrected[valid_corr])
            corrected[valid_corr] *= scale

        n_valid = int(valid.sum())
        n_ideal = ideal_idx.size
        n_nan = int(np.isnan(corrected).sum())
        ideal_pct = n_ideal / max(n_valid, 1) * 100
        nan_pct = n_nan / max(n, 1) * 100
        logging.info(
            f"    {stage_name:<5s}  {n_ideal:>8d}/{n_valid} ({ideal_pct:5.1f}%) fit,  "
            f"{n_nan:>8d}/{n} ({nan_pct:5.1f}%) NaN"
        )
        return corrected, rmse

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)
    grid_01 = np.linspace(0, 1, grid_size)
    extra = (mappability >= min_mappability) if mappability is not None else None
    cor, gc_rmse = _apply_stage(
        reads,
        gc,
        gc_lo,
        gc_hi,
        grid_01,
        seed=seed,
        stage_name="GC",
        extra_ideal_mask=extra,
    )

    if mappability is not None:
        map_lo = np.nanquantile(mappability, doutlier)
        map_hi = np.nanquantile(mappability, 1.0 - doutlier)
        cor, _ = _apply_stage(
            cor,
            mappability,
            map_lo,
            map_hi,
            grid_01,
            seed=seed + 1,
            stage_name="MAP",
        )

    if repliseq is not None:
        repli_finite = np.isfinite(repliseq)
        repli_lo = np.nanquantile(repliseq[repli_finite], doutlier)
        repli_hi = np.nanquantile(repliseq[repli_finite], 1.0 - doutlier)
        grid_repli = np.linspace(repli_lo, repli_hi, grid_size)
        cor, _ = _apply_stage(
            cor,
            repliseq,
            repli_lo,
            repli_hi,
            grid_repli,
            seed=seed + 2,
            stage_name="REPLI",
        )

    return cor.astype(np.float32), gc_rmse


def correct_readcount_quadreg(
    reads,
    gc,
    mappability=None,
    repliseq=None,
    doutlier=0.001,
    min_mappability=0.9,
    eps_quantile=0.01,
):
    """Quadratic median quantile regression bias correction.

    Fits RD ~ GC + GC**2 via median quantile regression on valid bins,
    then divides raw depth by predicted and rescales to preserve median.
    When ``repliseq`` is provided, RT + RT**2 terms are added to the model.
    Mappability is used as a filter (not a covariate).

    Returns
    -------
    np.ndarray
        Corrected read counts. NaN where correction is not possible.
    float
        RMSE from the GC fit.
    """
    import pandas as pd
    import statsmodels.formula.api as smf

    reads = reads.astype(np.float64)
    n = len(reads)

    gc_lo = np.nanquantile(gc, doutlier)
    gc_hi = np.nanquantile(gc, 1.0 - doutlier)

    valid = (reads > 0) & np.isfinite(gc) & (gc >= gc_lo) & (gc <= gc_hi)
    if mappability is not None:
        valid &= mappability >= min_mappability

    if repliseq is not None:
        valid &= np.isfinite(repliseq)

    if repliseq is not None:
        fit_df = pd.DataFrame(
            {"RD": reads[valid], "GC": gc[valid], "RT": repliseq[valid]}
        )
        pred_df = pd.DataFrame({"GC": gc, "RT": np.nan_to_num(repliseq, nan=0.0)})
        formula = "RD ~ GC + I(GC**2) + RT + I(RT**2)"
    else:
        fit_df = pd.DataFrame({"RD": reads[valid], "GC": gc[valid]})
        pred_df = pd.DataFrame({"GC": gc})
        formula = "RD ~ GC + I(GC**2)"

    n_fit = len(fit_df)
    logging.info(f"    MEDIAN  {n_fit:>8d}/{n} fitting bins, formula: {formula}")

    if n_fit < 10:
        logging.warning("    MEDIAN: too few fitting bins; skipping correction")
        return reads.astype(np.float32), 0.0

    res = smf.quantreg(formula, data=fit_df).fit(q=0.5)
    predicted = res.predict(pred_df).to_numpy()

    all_valid = (reads > 0) & np.isfinite(gc)
    n_all_valid = int(all_valid.sum())
    rmse = (
        float(np.sqrt(np.mean((reads[all_valid] - predicted[all_valid]) ** 2)))
        if n_all_valid > 0
        else 0.0
    )

    eps = (
        np.nanquantile(predicted[predicted > 0], eps_quantile)
        if (predicted > 0).any()
        else 1.0
    )
    den = np.clip(np.nan_to_num(predicted, nan=eps), eps, None)

    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = np.where(reads > 0, reads / den, np.nan)

    if mappability is not None:
        corrected[mappability < min_mappability] = np.nan

    valid_corr = np.isfinite(corrected) & (reads > 0)
    if valid_corr.any():
        scale = np.median(reads[valid_corr]) / np.median(corrected[valid_corr])
        corrected[valid_corr] *= scale

    n_nan = int(np.isnan(corrected).sum())
    logging.info(f"    MEDIAN  {n_nan:>8d}/{n} ({n_nan / max(n, 1) * 100:5.1f}%) NaN")

    return corrected.astype(np.float32), rmse
