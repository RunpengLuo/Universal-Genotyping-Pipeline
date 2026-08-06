"""Phasing helpers: phase-set grouping, phase application, flips, switch probabilities.

Everything downstream of the phaser that decides how a haplotype is carried across
SNPs and blocks: which SNPs share a phase set, how a per-SNP phase label turns REF/ALT
counts into A/B counts, where a phase flip splits a block, and the probability that a
haplotype switches between consecutive blocks.
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import issparse
from scipy.stats import beta as beta_dist


def setup_phaseset_groups(snps):
    """Ensure ``seg_id``/``PS`` columns exist; return grouping ``[region_id, seg_id, PS]``.

    ``seg_id`` (the breakpoint chunk from build_segment_bed) is the hard bin boundary;
    binning groups by it and never merges across it, while ``region_id`` (the arm) is
    carried for RDR/QC. When ``seg_id`` is absent (no global BED), it falls back to
    ``region_id`` so grouping is identical to the pre-seg_id behavior. If ``PS`` is
    absent, set ``PS=1``; when present every SNP must carry a non-null value.
    """
    assert "region_id" in snps.columns, "invalid SNP file"
    if "seg_id" not in snps.columns:
        snps["seg_id"] = snps["region_id"]
    if "PS" not in snps.columns:
        logging.info("PS not in SNP columns, setting PS=1 for all SNPs")
        snps["PS"] = 1
    else:
        assert snps["PS"].notna().all(), "unexpected SNP without PS in phased VCF"
    logging.info(
        f"#seg_id={snps['seg_id'].nunique()}, #phaseset={snps['PS'].nunique()}"
    )
    return ["region_id", "seg_id", "PS"]


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
    snps, a_mtx, b_mtx, grp_cols, tumor_sidx=0, epsilon=0.05, alpha=0.05
):
    """Detect phase flips between consecutive SNPs using Beta credible intervals.

    For each pair of consecutive SNPs within a group, compute a 95% Beta(b+1, a+1)
    credible interval for the BAF. If any tumor sample shows the two intervals
    confidently on opposite sides of a dead zone around 0.5, mark a phase boundary.

    Parameters
    ----------
    snps : pd.DataFrame
        SNP DataFrame with grouping columns.
    a_mtx : (N, M) ndarray
        A-allele counts (N SNPs, M samples).
    b_mtx : (N, M) ndarray
        B-allele counts (N SNPs, M samples).
    grp_cols : list of str
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
        Globally unique phase_group IDs aligned to snps index.
    """
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

    phase_group = np.zeros(len(snps), dtype=np.int64)
    global_pg = 0
    n_boundaries = 0
    n_groups_split = 0
    flip_records = []

    for _, grp in snps.groupby(grp_cols, sort=False):
        idx = grp.index.to_numpy()
        if len(idx) < 2:
            phase_group[idx] = global_pg
            global_pg += 1
            continue

        # Vectorized: check all consecutive pairs × all samples at once
        hi_prev, lo_curr = ci_hi[idx[:-1]], ci_lo[idx[1:]]
        lo_prev, hi_curr = ci_lo[idx[:-1]], ci_hi[idx[1:]]
        is_flip = (
            ((hi_prev < 0.5 - epsilon) & (lo_curr > 0.5 + epsilon))
            | ((lo_prev > 0.5 + epsilon) & (hi_curr < 0.5 - epsilon))
        ).any(axis=1)

        local_pg = np.concatenate([[0], np.cumsum(is_flip)])
        phase_group[idx] = global_pg + local_pg

        n_flip = int(is_flip.sum())
        n_boundaries += n_flip
        if n_flip > 0:
            n_groups_split += 1
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

        global_pg += int(local_pg[-1]) + 1

    logging.info(
        f"detect_phase_flips: epsilon={epsilon}, alpha={alpha}, "
        f"boundaries={n_boundaries}, groups_split={n_groups_split}, "
        f"new_phase_groups={global_pg}"
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

    return pd.Series(phase_group, index=snps.index, dtype=np.int64)


def interp_cM_blocks(
    blocks: pd.DataFrame,
    snp_info: pd.DataFrame,
    genetic_map: pd.DataFrame,
    block_id_col: str = "bbc_id",
):
    """Interpolate centimorgan distances between consecutive blocks using a genetic map.

    Parameters
    ----------
    blocks : pd.DataFrame
        Block-level DataFrame with a block ID column and ``#CHR`` column.
    snp_info : pd.DataFrame
        SNP DataFrame with a matching block ID column and ``POS`` column.
    genetic_map : pd.DataFrame
        Genetic map with ``#CHR``, ``POS``, and ``cM`` columns.
    block_id_col : str
        Name of the block ID column in both *blocks* and *snp_info*.

    Returns
    -------
    np.ndarray
        Inter-block cM distances (first block per chromosome gets 0).
    """
    blocks = blocks.copy(deep=True)
    blocks["dist_cM"] = 0.0

    hb_pos = snp_info.groupby(block_id_col, sort=False)["POS"].agg(
        snp_start="min", snp_end="max"
    )
    blocks = blocks.join(hb_pos, on=block_id_col)

    genetic_map_chrs = genetic_map.groupby(by="#CHR", sort=False, observed=True)
    for ch, ch_blocks in blocks.groupby(by="#CHR", sort=False, observed=True):
        ch_map = genetic_map_chrs.get_group(ch)
        start_cMs = np.interp(
            ch_blocks["snp_start"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )
        end_cMs = np.interp(
            ch_blocks["snp_end"].to_numpy(),
            ch_map["POS"].to_numpy(),
            ch_map["cM"].to_numpy(),
        )

        dist_cM = np.zeros(len(ch_blocks), dtype=np.float32)
        dist_cM[1:] = start_cMs[1:] - end_cMs[:-1]
        blocks.loc[ch_blocks.index, "dist_cM"] = np.maximum(dist_cM, 0.0)
    return blocks["dist_cM"].to_numpy()


def estimate_switchprobs_cM(dist_cms: np.ndarray, nu=1, min_switchprob=1e-6):
    """Convert cM distances to haplotype switch probabilities using the Haldane mapping function.

    Computes ``(1 - exp(-2 * nu * d)) / 2``, clipped at *min_switchprob*.

    Parameters
    ----------
    dist_cms : np.ndarray
        Inter-SNP or inter-block centimorgan distances.
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


def estimate_switchprobs_PS(blocks: pd.DataFrame, switchprob_ps=0.05):
    """Assign switch probabilities based on PS phaseset membership.

    Within the same phaseset, the probability is *switchprob_ps*; across
    different phasesets it is approximately 0.5.

    Parameters
    ----------
    blocks : pd.DataFrame
        DataFrame with a ``PS`` column indicating phaseset IDs.
    switchprob_ps : float
        Switch probability within the same phaseset.

    Returns
    -------
    np.ndarray
        Switch probabilities per block.
    """
    switch_bias = 1e-4
    same_block = blocks["PS"] == blocks["PS"].shift(1).fillna(False)
    switchprobs = np.where(
        same_block,
        switchprob_ps,  # within same PS phase block
        0.5 - switch_bias,  # across different PS phase block
    )
    return switchprobs
