"""Within-bb re-orientation of SNP phases by a naive-Bayes EM over the observations.

Last update: 2026-09-25

Functions:
- correct_bin_phases: per-SNP flip vector from a within-bb EM over the tumor observations
- fit_bin_phases / flips_from_fit: the same in two steps, so a gate sweep needs one fit
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.special import betaln, expit, logsumexp

from segmentation_utils import sum_observations_to_pseudobulk
from utils import log_hist


def _bb_logpmf(b, a, theta, tau):
    """Beta-binomial log-likelihood of *b* successes, up to a term constant in theta."""
    alpha, beta = theta * tau, (1.0 - theta) * tau
    return betaln(b + alpha, a + beta) - betaln(alpha, beta)


def _emission_table(b, t, thetas, tau):
    """(n_snps, n_obs, n_grid) float32 log-likelihood of the observed B count.

    Terms constant in *thetas* are dropped: the E-step uses a ratio between the two
    orientations and the M-step compares grid points, so the binomial coefficient
    cancels in both. *thetas* must be symmetric about 0.5, which lets the flipped
    orientation be read off the same table in reverse (``BetaBinom(t - b; t, theta) ==
    BetaBinom(b; t, 1 - theta)``) instead of stored twice.

    The cancellation inside :func:`_bb_logpmf` happens in float64 before the cast, so a
    large *tau* - where the two ``betaln`` terms are individually huge and nearly equal -
    still lands on the binomial limit to about nine digits.
    """
    n_snp, n_obs = b.shape
    out = np.empty((n_snp, n_obs, len(thetas)), dtype=np.float32)
    a = t - b
    for g, th in enumerate(thetas):
        out[:, :, g] = _bb_logpmf(b, a, th, tau)
    return out


def _bb_reduce(indptr, weighted):
    """Sum rows of *weighted* within each bb, given a sorted-by-bb row layout."""
    return np.add.reduceat(weighted, indptr[:-1], axis=0)


def _em_bins(b, t, bb_ids, tau, n_grid, grid_eps, n_restarts, max_iter, tol):
    """Run the naive-Bayes EM on every bb at once.

    All bbs share the iteration structure, so the E-step is one pass over the SNP axis
    and the M-step is one segment-sum per (observation, theta grid point). The emission
    table is built once and reused across restarts, which is what makes the restarts
    cheap. ``tau`` is fixed, not fitted: see :func:`correct_bin_phases`.

    Returns:
        ``(post, llr, n_bb)``: *post* is the per-SNP posterior that the SNP keeps its
        input orientation, *llr* the per-bb log-likelihood ratio against the balanced
        null.
    """
    n_grid = int(n_grid)
    thetas = np.linspace(float(grid_eps), 1.0 - float(grid_eps), n_grid)
    if n_grid % 2 == 0:
        raise ValueError("n_grid must be odd so the grid contains 0.5")
    rev = n_grid - 1 - np.arange(n_grid)
    mid = n_grid // 2

    order = np.argsort(bb_ids, kind="stable")
    b, t, bb_sorted = b[order], t[order], bb_ids[order]
    n_bb = int(bb_sorted.max()) + 1 if len(bb_sorted) else 0
    indptr = np.r_[0, np.flatnonzero(np.diff(bb_sorted)) + 1, len(bb_sorted)]
    seg_bb = bb_sorted[indptr[:-1]]
    seg_of_snp = np.searchsorted(seg_bb, bb_sorted)
    n_seg, n_obs = len(seg_bb), b.shape[1]

    tau = float(tau)
    logging.info(
        f"within-bb EM: {len(b)} SNPs, {n_seg} non-empty bbs of {n_bb}, {n_obs} "
        f"observations, tau={tau:g}, theta grid={n_grid}, restarts={n_restarts}"
    )

    table = _emission_table(b, t, thetas, tau)
    null_ll = _bb_reduce(indptr, table[:, :, mid].astype(np.float64)).sum(axis=1)
    best_ll = np.full(n_seg, -np.inf)
    best_post = np.zeros(len(b))
    # the null theta=0.5 is a fixed point of the EM and is always one of the restarts,
    # so the fit can never score below it and the reported llr stays non-negative
    for init in np.r_[np.linspace(0.05, 0.45, int(n_restarts)), 0.5]:
        k = np.full((n_seg, n_obs), int(np.argmin(np.abs(thetas - init))))
        prev = np.full(n_seg, -np.inf)
        for _ in range(int(max_iter)):
            ks = k[seg_of_snp]
            lk = np.take_along_axis(table, ks[:, :, None], axis=2)[:, :, 0].sum(axis=1)
            lf = np.take_along_axis(table, rev[ks][:, :, None], axis=2)[:, :, 0].sum(
                axis=1
            )
            post = expit(lk - lf)
            ll = _bb_reduce(
                indptr, logsumexp(np.c_[lk, lf], axis=1)[:, None] - np.log(2.0)
            ).ravel()
            w32 = post.astype(np.float32)
            keep_sum = _bb_reduce(indptr, w32[:, None, None] * table)
            flip_sum = _bb_reduce(indptr, (1.0 - w32)[:, None, None] * table)
            k = np.argmax(keep_sum + flip_sum[:, :, rev], axis=2)
            if np.all(ll - prev < float(tol)):
                prev = ll
                break
            prev = ll
        better = prev > best_ll
        best_ll = np.where(better, prev, best_ll)
        take = better[seg_of_snp]
        best_post[take] = post[take]

    llr = np.zeros(n_bb)
    llr[seg_bb] = best_ll - null_ll
    post_out = np.empty(len(b))
    post_out[order] = best_post
    return post_out, llr, n_bb


def fit_bin_phases(
    snps_bb: pd.DataFrame,
    a_mtx,
    b_mtx,
    tumor_cols,
    tau,
    *,
    n_grid,
    grid_eps,
    n_restarts,
    max_iter,
    tol,
    min_snps,
    obs_cluster_ids=None,
    n_obs_clusters=None,
) -> dict:
    """Run the within-bb EM and return its state, before the gate is applied.

    Split out of :func:`correct_bin_phases` because the gate is pure post-processing:
    one fit serves a whole ``min_llr`` sweep.

    Args:
        snps_bb: One row per SNP, carrying ``bb_id``; see :func:`correct_bin_phases`.
        a_mtx: (n_snps, n_observations) A-allele counts, dense or CSR.
        b_mtx: (n_snps, n_observations) B-allele counts, dense or CSR.
        tumor_cols: Columns the EM is fitted on.
        tau: Beta-binomial dispersion.
        n_grid: Points on the ``theta`` grid; must be odd so the grid contains 0.5.
        grid_eps: Distance of the grid's ends from 0 and 1.
        n_restarts: Initial ``theta`` values tried, besides the null at 0.5.
        max_iter: EM iterations per restart.
        tol: Per-bb log-likelihood gain below which a restart stops.
        min_snps: SNPs a bb needs before it is fitted at all.
        obs_cluster_ids: Single-cell only; pseudobulk grouping of the columns.
        n_obs_clusters: Number of those clusters.

    Every solver setting is required and carries no default here: each belongs to the
    run's configuration, so that ``config.yaml`` is the only place any of them is set.

    Returns:
        ``{"post", "usable", "bb_ids", "reads", "llr", "n_bb", "tau"}``. *post* is the
        per-SNP posterior of keeping the input orientation over the usable SNPs; *llr*
        is ``loglik - loglik at theta=0.5`` per bb, both marginal over the per-SNP
        orientation.
    """
    tumor_cols = np.asarray(tumor_cols, dtype=np.int64)
    a = a_mtx[:, tumor_cols]
    b = b_mtx[:, tumor_cols]
    if issparse(a):
        a, b = a.toarray(), b.toarray()
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    if obs_cluster_ids is not None:
        a = sum_observations_to_pseudobulk(a, obs_cluster_ids, n_obs_clusters)
        b = sum_observations_to_pseudobulk(b, obs_cluster_ids, n_obs_clusters)
    t = a + b

    bb_ids = snps_bb["bb_id"].to_numpy(np.int64)
    usable = (t.sum(axis=1) > 0) & (bb_ids >= 0)
    counts = np.bincount(bb_ids[usable], minlength=int(bb_ids.max()) + 1)
    usable &= counts[np.clip(bb_ids, 0, len(counts) - 1)] >= int(min_snps)
    if not usable.any():
        logging.warning("within-bb EM: no bb with enough covered SNPs")
        return dict(
            post=np.zeros(0),
            usable=usable,
            bb_ids=np.zeros(0, np.int64),
            reads=np.zeros(0),
            llr=np.zeros(0),
            n_bb=0,
        )

    post, llr, n_bb = _em_bins(
        b[usable],
        t[usable],
        bb_ids[usable],
        tau,
        n_grid,
        grid_eps,
        n_restarts,
        max_iter,
        tol,
    )
    return dict(
        post=post,
        usable=usable,
        bb_ids=bb_ids[usable],
        reads=t[usable].sum(axis=1),
        llr=llr,
        n_bb=n_bb,
        tau=float(tau),
    )


def flips_from_fit(fit: dict, n_snps: int, min_llr) -> np.ndarray:
    """Turn a :func:`fit_bin_phases` state into a per-SNP flip vector.

    Flipping every SNP in a bb is an exact symmetry of the likelihood, so the EM leaves
    each bb's global orientation free. It is resolved towards the input phasing: of the
    two equivalent answers, the one that re-orients at most half of the bb's reads. The
    phaser's cross-bb frame therefore survives and the bb-level ``switchprobs`` keep
    their meaning. A run with no phaser has an arbitrary but deterministic REF/ALT frame
    and anchoring to it is equally harmless, since its ``switchprobs`` are ~0.5 anyway.

    Args:
        fit: The dict :func:`fit_bin_phases` returned.
        n_snps: Length of the output, i.e. the number of rows in the SNP frame.
        min_llr: Evidence a bb must show of having any allelic imbalance before its
            fitted orientation is accepted, as the log of a likelihood ratio; see
            :func:`correct_bin_phases`. ``0`` accepts every bb.

    Returns:
        (n_snps,) int8, 1 where the SNP's A and B counts must be swapped.
    """
    flips = np.zeros(n_snps, dtype=np.int8)
    if fit["n_bb"] == 0:
        return flips
    post, sub_ids, reads = fit["post"], fit["bb_ids"], fit["reads"]
    n_bb = fit["n_bb"]
    hard = post < 0.5
    tsum = np.bincount(sub_ids, weights=reads, minlength=n_bb)
    moved = np.bincount(sub_ids, weights=hard * reads, minlength=n_bb)
    global_flip = moved > tsum - moved

    gated = fit["llr"] < float(min_llr)
    sub_flip = (hard ^ global_flip[sub_ids]) & ~gated[sub_ids]
    flips[np.flatnonzero(fit["usable"])] = sub_flip.astype(np.int8)
    uniq = np.unique(sub_ids)
    log_hist(fit["llr"][uniq], "per-bb evidence of imbalance (log likelihood ratio)")
    logging.info(
        f"within-bb EM: min_llr={min_llr}, "
        f"{int((~gated[uniq]).sum())}/{len(uniq)} bbs above the gate, "
        f"{int(flips.sum())}/{n_snps} SNPs re-oriented"
    )
    return flips


def correct_bin_phases(
    snps_bb: pd.DataFrame,
    a_mtx,
    b_mtx,
    tumor_cols,
    tau,
    min_llr,
    *,
    n_grid,
    grid_eps,
    n_restarts,
    max_iter,
    tol,
    min_snps,
    obs_cluster_ids=None,
    n_obs_clusters=None,
) -> np.ndarray:
    """Per-SNP flip vector re-orienting SNP phases inside each bb.

    The model is a naive Bayes EM over the tumor observations, after HATCHet2's
    ``multisample_em`` (``hatchet/utils/combine_counts.py``,
    https://github.com/raphael-group/hatchet, doi:10.1038/s41587-020-0661-6); Alleloscope
    (doi:10.1038/s41587-021-00911-w) uses the same shape of per-region latent-phase EM.
    The per-SNP phase latent is an independent Bernoulli - there is no transition prior,
    no genetic map and no switch probability here, because the coupling that identifies
    the orientation runs across *observations*, not across position. A SNP whose
    orientation is ambiguous in one tumor is resolved by the others.

    Why this exists: ``apply_phase_to_mat`` orients A/B from one per-SNP ``PHASE`` bit
    and nothing re-orients within a bb, so a bb spanning a phaser switch error sums two
    anti-phased halves and averages its BAF back toward 0.5.

    Emission: ``BetaBinom(b; t, tau*theta, tau*(1-theta))`` per tumor observation, where
    ``theta`` is that bb's BAF in that observation. A binomial (``tau -> inf``) reads the
    overdispersion of real data - reference and capture bias - as evidence of a flip.
    ``tau`` is ONE number for the whole genome and is **fixed, not fitted**: it is not
    identifiable jointly with ``theta``, which is free per bb and per observation and so
    already absorbs the variation between bbs, leaving the residual scatter of a SNP
    around its own bb's ``theta`` binomial. Smaller ``tau`` re-orients fewer SNPs.

    ``theta`` is maximized over a grid rather than in closed form; the grid is symmetric
    about 0.5, so the flipped orientation is the same table read backwards.

    Each bb's global orientation is resolved towards the input phasing rather than folded
    to BAF <= 0.5, so the phaser's cross-bb frame survives; see :func:`flips_from_fit`.

    Gating. A balanced bb has no orientation to find, but the EM will still split its
    SNPs and report an imbalance - the classic latent-variable overfit. It inflates the
    spread of BAF around 0.5 far more than it shifts the median, because the EM picks the
    partition that *maximizes* imbalance, so a minority of bbs land a long way out.

    Every bb is fitted; ``min_llr`` then decides whether that fit is kept. It is the
    evidence a bb must show of being imbalanced at all:
    ``loglik - loglik at theta=0.5``. Both sides marginalize over the per-SNP
    orientation under a uniform prior, so they are likelihoods of the same data under two
    nested models that differ only in the free BAF per observation. At the null that
    marginalization is a no-op - the two orientations of a SNP have identical emissions
    when ``theta = 0.5`` - which is also why the null makes orientation unidentifiable. A bb below the floor keeps its
    input orientation untouched. It is a log likelihood ratio in nats, so a floor of
    ``L`` demands the imbalanced model be ``exp(L)`` times as likely as the balanced one.
    It is NOT a chi-square statistic - that would be twice this, and in any case the
    alternative carries one free BAF per observation *and* a latent orientation per SNP,
    so Wilks does not apply and the floor has to be calibrated rather than read off a
    table. Because evidence accumulates with reads, one floor is automatically stricter
    for a small bb, which is where the overfit lives. Calibrate on regions known to be
    balanced: raise it until their BAF spread returns to what it was before the EM ran.

    Args:
        snps_bb: One row per SNP, carrying ``bb_id``. Rows are aligned to the matrix rows
            but are not required to be sorted; this function groups by ``bb_id`` itself
            and never permutes the matrices.
        a_mtx: (n_snps, n_observations) A-allele counts, dense or CSR. Never mutated.
        b_mtx: (n_snps, n_observations) B-allele counts, dense or CSR. Never mutated.
        tumor_cols: Columns the EM is fitted on. A matched normal is excluded: its true
            BAF is 0.5, so it carries no orientation signal. No tumor column means no
            SNP is re-oriented.
        obs_cluster_ids: Single-cell only; each column's observation cluster, so per-cell
            columns are pseudobulked before the fit. Per-cell counts are far too sparse.
        tau: Beta-binomial dispersion; see above.
        min_llr: Evidence floor for accepting a bb's fitted orientation; see above.
        n_grid, grid_eps, n_restarts, max_iter, tol, min_snps: Solver settings, passed
            through to :func:`fit_bin_phases`.
        n_obs_clusters: Number of those clusters.

    No argument carries a default: every one of them is a run configuration key, so
    ``config.yaml`` is the only place any of them is set.

    Returns:
        (n_snps,) int8, 1 where the SNP's A and B counts must be swapped, positionally
        aligned to *snps_bb* and to the matrix rows.
    """
    n_snp = len(snps_bb)
    if len(np.asarray(tumor_cols)) == 0:
        logging.warning(
            "within-bb EM: no tumor observation to fit on, no SNP re-oriented"
        )
        return np.zeros(n_snp, dtype=np.int8)

    fit = fit_bin_phases(
        snps_bb,
        a_mtx,
        b_mtx,
        tumor_cols,
        tau,
        n_grid=n_grid,
        grid_eps=grid_eps,
        n_restarts=n_restarts,
        max_iter=max_iter,
        tol=tol,
        min_snps=min_snps,
        obs_cluster_ids=obs_cluster_ids,
        n_obs_clusters=n_obs_clusters,
    )
    return flips_from_fit(fit, n_snp, min_llr)
