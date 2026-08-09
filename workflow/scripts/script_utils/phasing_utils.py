"""Phasing helpers: phase clusters, phase application, flips, switch probabilities.

Everything downstream of the phaser that decides how a haplotype is carried across SNPs
and bbs: which SNPs share a phase cluster (the VCF ``PS`` tag), how a per-SNP phase
label turns REF/ALT counts into A/B counts, where a phase flip splits a cluster, and the
probability that a haplotype switches between consecutive bbs.
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import beta as beta_dist


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


def detect_phase_flips(
    snps, a_mtx, b_mtx, cluster_cols, tumor_sidx=0, epsilon=0.05, alpha=0.05
):
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
    a_mtx : (N, M) ndarray
        A-allele counts (N SNPs, M samples).
    b_mtx : (N, M) ndarray
        B-allele counts (N SNPs, M samples).
    cluster_cols : list of str
        Columns to group SNPs by (e.g. ["region_id", "PS"]).
    tumor_sidx : int
        Index of first tumor sample column.
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

    a_tumor = (
        a_mtx[:, tumor_sidx:].toarray() if issparse(a_mtx) else a_mtx[:, tumor_sidx:]
    ).astype(np.float64)
    b_tumor = (
        b_mtx[:, tumor_sidx:].toarray() if issparse(b_mtx) else b_mtx[:, tumor_sidx:]
    ).astype(np.float64)

    ci_lo = beta_dist.ppf(alpha / 2, b_tumor + 1, a_tumor + 1)
    ci_hi = beta_dist.ppf(1 - alpha / 2, b_tumor + 1, a_tumor + 1)

    # Observed BAF for logging
    tot_tumor = a_tumor + b_tumor
    baf = np.divide(
        b_tumor, tot_tumor, out=np.full_like(b_tumor, np.nan), where=tot_tumor > 0
    )

    phase_cluster = np.zeros(len(snps), dtype=np.int64)
    global_pc = 0
    n_boundaries = 0
    n_clusters_split = 0
    flip_records = []

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

        n_flip = int(is_flip.sum())
        n_boundaries += n_flip
        if n_flip > 0:
            n_clusters_split += 1
            for fp in np.where(is_flip)[0]:
                pi, ci = idx[fp], idx[fp + 1]
                mean_diff = np.nanmean(np.abs(baf[pi] - baf[ci]))
                flip_records.append(
                    (
                        snps.iat[pi, snps.columns.get_loc("#CHR")],
                        snps.iat[pi, snps.columns.get_loc("POS0")],
                        snps.iat[ci, snps.columns.get_loc("POS0")],
                        baf[pi],
                        baf[ci],
                        mean_diff,
                    )
                )

        global_pc += int(local_pc[-1]) + 1

    logging.info(
        f"detect_phase_flips: epsilon={epsilon}, alpha={alpha}, "
        f"boundaries={n_boundaries}, clusters_split={n_clusters_split}, "
        f"new_phase_groups={global_pc}"
    )

    if flip_records:
        flip_records.sort(key=lambda r: r[5], reverse=True)
        n_show = min(10, len(flip_records))
        logging.info(f"top {n_show} flips by |ΔBAF| (largest gap):")
        for chrom, p1, p2, b1, b2, d in flip_records[:n_show]:
            b1s = ",".join(f"{v:.3f}" for v in b1)
            b2s = ",".join(f"{v:.3f}" for v in b2)
            logging.info(f"  {chrom}:{p1}-{p2}  BAF=[{b1s}]→[{b2s}]  |Δ|={d:.3f}")
        if len(flip_records) > n_show:
            logging.info(f"bottom {n_show} flips by |ΔBAF| (smallest gap):")
            for chrom, p1, p2, b1, b2, d in flip_records[-n_show:]:
                b1s = ",".join(f"{v:.3f}" for v in b1)
                b2s = ",".join(f"{v:.3f}" for v in b2)
                logging.info(f"  {chrom}:{p1}-{p2}  BAF=[{b1s}]→[{b2s}]  |Δ|={d:.3f}")

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


def estimate_switchprobs_PS(bbs: pd.DataFrame, switchprob_ps=0.05):
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
    switch_bias = 1e-4
    same_ps = bbs["PS"] == bbs["PS"].shift(1).fillna(False)
    switchprobs = np.where(
        same_ps,
        switchprob_ps,  # within the same phase cluster
        0.5 - switch_bias,  # across phase clusters
    )
    return switchprobs
