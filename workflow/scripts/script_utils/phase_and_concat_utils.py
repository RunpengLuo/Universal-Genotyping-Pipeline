"""SNP and allele helpers for the phase_and_concat scripts.

Three groups, in pipeline order:

1. allele matrices - per-replicate counts onto the parent SNP set, then stacked
2. SNP ranges      - the per-SNP ``[START, END)`` split over the regions
3. SNP filters     - one bool mask each, combined by ``apply_masks_to_df``

The downstream half of the pipeline is ``combine_counts_utils``; nothing is shared.
"""

import logging

import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix, hstack
from scipy.stats import beta

from io_utils import read_BED
from range_utils import assign_pos_to_range, overlaps_any_range


##################################################
# allele matrices: per-replicate counts -> one joint matrix


def map_allele_mat_to_snps(
    parent_keys: pd.Index,
    child_snps: pd.DataFrame,
    tot_mtx: csr_matrix,
    ad_mtx: csr_matrix,
    ncells: int,
):
    """Map a replicate's DP/AD matrices onto the parent SNP set (0-fill missing loci).

    Args:
        parent_keys: parent SNP KEYs (``#CHROM_POS``) in output row order.
        child_snps: DataFrame with ``KEY`` + ``RAW_SNP_DF_IDX`` for this replicate's loci
            (in the matrices' row order).
        tot_mtx, ad_mtx: ``(m, ncells)`` count matrices, ``m == len(child_snps)``.
        ncells: number of observations (1 for bulk pseudobulk; #barcodes for single-cell).

    Returns:
        ``(tot_canon, ad_canon)`` csr matrices of shape ``(len(parent_keys), ncells)``.
    """
    M = len(parent_keys)
    m = len(child_snps)

    tot_mtx = tot_mtx.tocsr()
    ad_mtx = ad_mtx.tocsr()
    assert tot_mtx.shape == ad_mtx.shape, (
        f"count matrices, shape mismatch: {tot_mtx.shape} vs {ad_mtx.shape}"
    )
    assert tot_mtx.shape == (m, ncells), (
        f"count matrix, shape {tot_mtx.shape}, expected {(m, ncells)}"
    )

    raw_snp_df_idx = child_snps["RAW_SNP_DF_IDX"].to_numpy()
    tot_mtx = tot_mtx[raw_snp_df_idx, :]
    ad_mtx = ad_mtx[raw_snp_df_idx, :]

    dup_mask = child_snps["KEY"].duplicated(keep=False).to_numpy()
    n_dup_snps = int(dup_mask.sum())
    if n_dup_snps > 0:
        n_dup_keys = int(child_snps.loc[dup_mask, "KEY"].nunique())
        logging.warning(
            f"#found duplicated SNP #CHR/POS, {n_dup_snps}/{m} SNPs over {n_dup_keys} keys, drop all"
        )
        child_snps = child_snps.loc[~dup_mask].reset_index(drop=True)
        tot_mtx = tot_mtx[~dup_mask, :]
        ad_mtx = ad_mtx[~dup_mask, :]

    m = len(child_snps)
    child_keys = pd.Index(child_snps["KEY"])

    child_loc = child_keys.get_indexer(parent_keys)
    present = child_loc >= 0
    logging.info(f"located SNPs in parent={present.sum()}/{M}")

    common_pidx = np.flatnonzero(present)
    common_cidx = child_loc[present]

    tot_present = tot_mtx[common_cidx, :].tocoo()
    tot_canon = csr_matrix(
        (tot_present.data, (common_pidx[tot_present.row], tot_present.col)),
        shape=(M, ncells),
        dtype=tot_mtx.dtype,
    )

    ad_present = ad_mtx[common_cidx, :].tocoo()
    ad_canon = csr_matrix(
        (ad_present.data, (common_pidx[ad_present.row], ad_present.col)),
        shape=(M, ncells),
        dtype=ad_mtx.dtype,
    )

    return tot_canon, ad_canon


def hstack_replicate_mats(tot_list: list, ad_list: list):
    """Horizontally stack per-replicate total and alt-count matrices.

    Derives the REF matrix as ``TOT - ALT``.

    Parameters
    ----------
    tot_list : list of csr_matrix
        Per-replicate total depth matrices (SNP features x observations).
    ad_list : list of csr_matrix
        Per-replicate alt-allele count matrices (SNP features x observations).

    Returns
    -------
    tuple[csr_matrix, csr_matrix, csr_matrix]
        ``(tot_mtx, ref_mtx, alt_mtx)`` concatenated across replicates.
    """
    if len(tot_list) == 1:
        return tot_list[0], tot_list[0] - ad_list[0], ad_list[0]
    tot_mtx = hstack(tot_list, format="csr")
    alt_mtx = hstack(ad_list, format="csr")
    ref_mtx = tot_mtx - alt_mtx
    return tot_mtx, ref_mtx, alt_mtx


##################################################
# SNP ranges: the per-SNP range split


def build_pos_ranges(snps: pd.DataFrame, regions: pd.DataFrame, colname="region_id"):
    """Split each region into one ``[START, END)`` range per position it contains.

    A position owns the span between the midpoints of its neighbours, bounded by the
    region edges; a lone position owns the whole region. Positions outside every region
    keep ``START = END = 0`` and an empty ``region_id``. Membership is by 0-based
    ``POS0``, like every other range operation.

    Args:
        snps: Positions with ``#CHR``, ``POS0``, sorted genomically.
        regions: Ranges with ``#CHR``, ``START``, ``END``, ``region_id`` and
            optionally ``seg_id``.
        colname: Column to write the region identifier into.

    Returns:
        *snps* with ``START``, ``END``, ``BLOCKSIZE``, *colname* (and ``seg_id``).
    """
    has_seg = "seg_id" in regions.columns
    regions = regions.reset_index(drop=True)
    regions["_reg_id"] = np.arange(len(regions))

    snps, _ = assign_pos_to_range(snps, regions, ref_id="_reg_id", pos_col="POS0")
    snps["START"] = 0
    snps["END"] = 0
    snps[colname] = ""
    if has_seg:
        snps["seg_id"] = ""

    inside = snps["_reg_id"].notna()
    n_out = int((~inside).sum())
    if n_out:
        logging.info(f"SNP ranges: {n_out}/{len(snps)} SNPs outside every region")

    reg_start = regions["START"].to_numpy()
    reg_end = regions["END"].to_numpy()
    for reg_id, grp in snps.loc[inside].groupby("_reg_id", sort=False):
        idx = grp.index.to_numpy()
        r = int(reg_id)
        snps.loc[idx, colname] = regions.at[r, "region_id"]
        if has_seg:
            snps.loc[idx, "seg_id"] = regions.at[r, "seg_id"]
        pos0 = grp["POS0"].to_numpy()
        if len(idx) == 1:
            bounds = np.array([reg_start[r], reg_end[r]])
        else:
            mids = np.ceil((pos0[:-1] + pos0[1:]) / 2).astype(np.uint32)
            bounds = np.concatenate([[reg_start[r]], mids, [reg_end[r]]])
        snps.loc[idx, "START"] = bounds[:-1]
        snps.loc[idx, "END"] = bounds[1:]

    snps.drop(columns="_reg_id", inplace=True)
    snps["BLOCKSIZE"] = snps["END"] - snps["START"]
    return snps


##################################################
# SNP filters: one bool mask each, combined by apply_masks_to_df


def apply_masks_to_df(df: pd.DataFrame, *masks):
    """Keep the rows of *df* that every mask keeps, reindexed from 0.

    Args:
        df: Frame to filter. Copied, not modified.
        *masks: One or more bool arrays of length ``len(df)``.

    Returns:
        ``(df, keep)`` - the surviving rows, ``reset_index(drop=True)``, and the
        combined mask over the INPUT rows, for subsetting a parallel matrix the
        same way.
    """
    assert masks, "apply_masks_to_df, no mask given"
    keep = np.logical_and.reduce([np.asarray(m, dtype=bool) for m in masks])
    assert len(keep) == len(df), f"mask length {len(keep)} != {len(df)} rows"
    return df.loc[keep].reset_index(drop=True), keep


def get_mask_by_region(snps: pd.DataFrame, regions: pd.DataFrame):
    """Keep SNPs inside any region."""
    mask = overlaps_any_range(snps, regions)
    logging.info(f"filter by region, #passed SNPs={np.sum(mask)}/{len(snps)}")
    return mask


def get_mask_by_blacklist(snps: pd.DataFrame, blacklist_bed):
    """Drop SNPs inside any blacklisted range; all-True when no blacklist is given."""
    if blacklist_bed is None:
        return np.ones(len(snps), dtype=bool)
    hit = overlaps_any_range(snps, read_BED(blacklist_bed))
    num_passes = len(snps) - np.sum(hit)
    logging.info(f"filter by blacklist, #passed SNPs={num_passes}/{len(snps)}")
    return ~hit


def get_mask_by_exon(snps: pd.DataFrame):
    """Keep only exonic SNPs.

    Reads the ``feature_type`` column that ``feature_utils.annotate_feature_type`` writes.
    The caller gates the filter on the ``exon_only`` config key.
    """
    is_exon = (snps["feature_type"] == "exon").to_numpy()
    logging.info(f"filter by exon, #passed SNPs={np.sum(is_exon)}/{len(snps)}")
    return is_exon


def get_mask_by_depth(snps: pd.DataFrame, tot_mtx: csr_matrix, min_dp=1):
    """Return a boolean mask keeping SNPs where every sample meets the minimum depth.

    Parameters
    ----------
    snps : pd.DataFrame
        SNP info DataFrame (used only for logging).
    tot_mtx : csr_matrix
        Total depth matrix (SNPs x samples).
    min_dp : int
        Minimum depth threshold per sample.

    Returns
    -------
    np.ndarray
        Boolean mask of length ``len(snps)``.
    """
    mask = np.all(tot_mtx >= min_dp, axis=1)
    logging.info(
        f"filter by depth, min_dp={min_dp}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


def get_mask_by_het_balanced(
    snps: pd.DataFrame,
    ref_mtx: csr_matrix,
    alt_mtx: csr_matrix,
    gamma: float,
    normal_idx=0,
):
    """
    mask SNPs if normal sample failed beta-posterior credible interval test with beta(1, 1) prior.
    """
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    q = np.array([p_lower, p_upper])
    het_cred_ints = beta.ppf(
        q[None, :],
        ref_mtx[:, normal_idx][:, None] + 1,
        alt_mtx[:, normal_idx][:, None] + 1,
    )
    mask = (het_cred_ints[:, 0] <= 0.5) & (0.5 <= het_cred_ints[:, 1])
    logging.info(
        f"filter by balanced Het-SNPs on normal sample, gamma={gamma}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask
