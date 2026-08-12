"""Carry a haplotype across SNPs and bbs, after the phaser.

Last update: 2026-08-11

Functions:
- apply_phase_to_mat: turn REF/ALT counts into phased A/B counts
- detect_phase_flips: split a cluster where the haplotype orientation switches
- interp_cM_between_bbs: centimorgan distance between consecutive bbs
- estimate_switchprobs_cM: Haldane switch probability from a cM distance
- estimate_switchprobs_PS: switch probability from phase-cluster (PS) membership
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import beta as beta_dist

from utils import log_hist


def apply_phase_to_mat(tot_mtx, ref_mtx, alt_mtx, phases):
    """Apply per-SNP phase labels to produce phased A/B allele count matrices.

    Parameters
    ----------
    tot_mtx : sparse or ndarray
        Total depth matrix (SNPs x cells/samples).
    ref_mtx : sparse or ndarray
        Reference allele count matrix.
    alt_mtx : sparse or ndarray
        Alternate allele count matrix.
    phases : np.ndarray
        Per-SNP phase labels (0 or 1).

    Returns
    -------
    tuple
        ``(a_mtx, b_mtx)`` — phased allele count matrices.
    """
    p = phases[:, None]
    if issparse(ref_mtx):
        b_mtx = ref_mtx.multiply(p) + alt_mtx.multiply(1 - p)
        b_mtx.data = np.rint(b_mtx.data).astype(np.int32)
    else:
        b_mtx = ref_mtx * p + alt_mtx * (1 - p)
        b_mtx = np.round(b_mtx).astype(np.int32)
    a_mtx = tot_mtx - b_mtx
    return a_mtx, b_mtx


def _baf(b, a):
    """B-allele fraction, NaN where the two counts sum to zero."""
    tot = a + b
    return np.divide(b, tot, out=np.full_like(b, np.nan), where=tot > 0)


def detect_phase_flips(snps, a_mtx, b_mtx, cluster_cols, epsilon=0.05, alpha=0.05):
    """Detect phase flips between neighbouring SNPs using Beta credible intervals.

    For each pair of genomically adjacent SNPs within a group, compute a 95%
    Beta(b+1, a+1) credible interval for the BAF. If any tumor sample shows the two
    intervals confidently on opposite sides of a dead zone around 0.5, mark a phase
    boundary. Each group is ordered by ``POS0`` here, so *snps* need not arrive sorted
    and the matrices are never permuted.

    Parameters
    ----------
    snps : pd.DataFrame
        SNP DataFrame with grouping columns.
    a_mtx, b_mtx : (N, M) dense ndarray
        A- and B-allele counts over the M TUMOR observations; the caller slices.
    cluster_cols : list of str
        Columns to group SNPs by (e.g. ["region_id", "PS"]).
    epsilon : float
        Half-width of dead zone around 0.5. Default 0.05 → dead zone [0.45, 0.55].
    alpha : float
        Significance level for credible intervals. Default 0.05 → 95% CI.

    Returns
    -------
    pd.Series
        Globally unique phase-cluster IDs aligned to the snps index.
    """
    orig_index = snps.index
    snps = snps.reset_index(drop=True)  # idx below indexes the matrices positionally
    pos0 = snps["POS0"].to_numpy()

    a_tumor = np.asarray(a_mtx, dtype=np.float64)
    b_tumor = np.asarray(b_mtx, dtype=np.float64)

    ci_lo = beta_dist.ppf(alpha / 2, b_tumor + 1, a_tumor + 1)
    ci_hi = beta_dist.ppf(1 - alpha / 2, b_tumor + 1, a_tumor + 1)

    phase_cluster = np.zeros(len(snps), dtype=np.int64)
    global_pc = 0
    n_clusters_split = 0
    flip_gaps = []

    for _, grp in snps.groupby(cluster_cols, sort=False):
        idx = grp.index.to_numpy()
        if len(idx) < 2:
            phase_cluster[idx] = global_pc
            global_pc += 1
            continue
        # a flip is between genomic neighbours, not between adjacent rows
        idx = idx[np.argsort(pos0[idx], kind="stable")]

        # Vectorized: check all consecutive pairs × all samples at once
        hi_prev, lo_curr = ci_hi[idx[:-1]], ci_lo[idx[1:]]
        lo_prev, hi_curr = ci_lo[idx[:-1]], ci_hi[idx[1:]]
        is_flip = (
            ((hi_prev < 0.5 - epsilon) & (lo_curr > 0.5 + epsilon))
            | ((lo_prev > 0.5 + epsilon) & (hi_curr < 0.5 - epsilon))
        ).any(axis=1)

        local_pc = np.concatenate([[0], np.cumsum(is_flip)])
        phase_cluster[idx] = global_pc + local_pc

        flips = np.flatnonzero(is_flip)
        if len(flips):
            n_clusters_split += 1
            prev, curr = idx[flips], idx[flips + 1]
            gap = np.abs(
                _baf(b_tumor[prev], a_tumor[prev]) - _baf(b_tumor[curr], a_tumor[curr])
            )
            flip_gaps.append(np.nanmean(gap, axis=1))

        global_pc += int(local_pc[-1]) + 1

    gaps = np.concatenate(flip_gaps) if flip_gaps else np.zeros(0)
    logging.info(
        f"detect_phase_flips: epsilon={epsilon}, alpha={alpha}, "
        f"boundaries={len(gaps)}, clusters_split={n_clusters_split}, "
        f"new_phase_groups={global_pc}"
    )
    log_hist(gaps, "|dBAF| across a detected flip")

    return pd.Series(phase_cluster, index=orig_index, dtype=np.int64)


def interp_cM_between_bbs(
    bbs: pd.DataFrame,
    snp_info: pd.DataFrame,
    genetic_map: pd.DataFrame,
    bb_id_col: str = "bb_id",
):
    """Interpolate centimorgan distances between consecutive bbs using a genetic map.

    Parameters
    ----------
    bbs : pd.DataFrame
        bb-level DataFrame with a bb ID column and ``#CHR`` column.
    snp_info : pd.DataFrame
        SNP DataFrame with a matching bb ID column and ``POS`` column.
    genetic_map : pd.DataFrame
        Genetic map with ``#CHR``, ``POS``, and ``cM`` columns.
    bb_id_col : str
        Name of the bb ID column in both *bbs* and *snp_info*.

    Returns
    -------
    np.ndarray
        Inter-bb cM distances (first bb per chromosome gets 0).
    """
    logging.info("interpolate phase switchprobs for bbs.")
    bbs = bbs.copy(deep=True)
    bbs["dist_cM"] = 0.0

    hb_pos = snp_info.groupby(bb_id_col, sort=False)["POS"].agg(
        snp_start="min", snp_end="max"
    )
    bbs = bbs.join(hb_pos, on=bb_id_col)
    # a bb holding no SNP falls back to its own span; a NaN here would interpolate to NaN
    # and propagate to the next bb through dist_cM
    bbs["snp_start"] = bbs["snp_start"].fillna(bbs["START"] + 1)
    bbs["snp_end"] = bbs["snp_end"].fillna(bbs["END"])

    genetic_map_chrs = genetic_map.groupby(by="#CHR", sort=False, observed=True)
    for ch, ch_bbs in bbs.groupby(by="#CHR", sort=False, observed=True):
        ch_map = genetic_map_chrs.get_group(ch)
        start_cMs = np.interp(
            ch_bbs["snp_start"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )
        end_cMs = np.interp(
            ch_bbs["snp_end"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )

        dist_cM = np.zeros(len(ch_bbs), dtype=np.float32)
        dist_cM[1:] = start_cMs[1:] - end_cMs[:-1]
        bbs.loc[ch_bbs.index, "dist_cM"] = np.maximum(dist_cM, 0.0)
    return bbs["dist_cM"].to_numpy()


def estimate_switchprobs_cM(dist_cms: np.ndarray, nu=1, min_switchprob=1e-6):
    """Convert cM distances to haplotype switch probabilities using the Haldane mapping function.

    Computes ``(1 - exp(-2 * nu * d)) / 2``, clipped at *min_switchprob*.

    Parameters
    ----------
    dist_cms : np.ndarray
        Inter-SNP or inter-bb centimorgan distances.
    nu : float
        Scaling factor for the Haldane function.
    min_switchprob : float
        Minimum switch probability (floor).

    Returns
    -------
    np.ndarray
        Switch probabilities.
    """
    switchprobs = (1 - np.exp(-2 * nu * dist_cms)) / 2.0
    return np.clip(switchprobs, a_min=min_switchprob, a_max=None)


def estimate_switchprobs_PS(bbs: pd.DataFrame, switchprob_ps=0.05, switch_bias=1e-4):
    """Assign switch probabilities based on phase-cluster (``PS``) membership.

    Within the same phase cluster the probability is *switchprob_ps*; across clusters
    it is approximately 0.5.

    Parameters
    ----------
    bbs : pd.DataFrame
        DataFrame with a ``PS`` column indicating the phase cluster.
    switchprob_ps : float
        Switch probability within the same phase cluster.

    Returns
    -------
    np.ndarray
        Switch probabilities per bb.
    """
    logging.info("assign phase switchprobs for bbs given PS")
    same_ps = bbs["PS"] == bbs["PS"].shift(1).fillna(False)
    switchprobs = np.where(
        same_ps,
        switchprob_ps,  # within the same phase cluster
        0.5 - switch_bias,  # across phase clusters
    )
    return switchprobs
