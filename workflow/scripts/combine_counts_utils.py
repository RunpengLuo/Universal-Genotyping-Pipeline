"""Utility functions for SNP/allele processing, phasing, and filtering.

Used by phase_and_concat.py, combine_counts.py, and combine_counts_nonbulk.py.
"""

import logging

import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta

from io_utils import read_VCF
from utils import stamp_path


def canon_mat_one_replicate(
    parent_keys: pd.Index,
    vcf_file: str,
    tot_mtx_file: str,
    ad_mtx_file: str,
    ncells: int,
):
    """
    Map SNP by barcode DP/AD sparse mats to same SNP position index defined by parent SNP file.
    """
    M = len(parent_keys)
    child_snps = read_VCF(vcf_file, addkey=True)
    m = len(child_snps)

    tot_mtx: csr_matrix = mmread(tot_mtx_file).tocsr()
    ad_mtx: csr_matrix = mmread(ad_mtx_file).tocsr()
    assert tot_mtx.shape == ad_mtx.shape
    assert tot_mtx.shape == (m, ncells)

    raw_snp_ids = child_snps["RAW_SNP_IDX"].to_numpy()
    tot_mtx = tot_mtx[raw_snp_ids, :]
    ad_mtx = ad_mtx[raw_snp_ids, :]

    dup_mask = child_snps["KEY"].duplicated(keep=False).to_numpy()
    n_dup_rows = int(dup_mask.sum())
    if n_dup_rows > 0:
        n_dup_keys = int(child_snps.loc[dup_mask, "KEY"].nunique())
        logging.warning(
            f"#found duplicated SNP #CHR/POS, #rows={n_dup_rows}/{m} #SNPs={n_dup_keys}, drop all"
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


def merge_mats(tot_list: list, ad_list: list):
    """Horizontally stack per-replicate total and alt-count matrices.

    Derives the REF matrix as ``TOT - ALT``.

    Parameters
    ----------
    tot_list : list of csr_matrix
        Per-replicate total depth matrices (SNPs x cells).
    ad_list : list of csr_matrix
        Per-replicate alt-allele count matrices (SNPs x cells).

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


def compute_af_per_sample(tot_mtx, b_mtx, i: int):
    """Compute per-SNP allele frequency for a single sample column.

    Parameters
    ----------
    tot_mtx : sparse or ndarray
        Total depth matrix (SNPs x samples).
    b_mtx : sparse or ndarray
        B-allele count matrix.
    i : int
        Sample column index.

    Returns
    -------
    np.ndarray
        Allele frequency per SNP; ``NaN`` where depth is zero.
    """
    tot_col = tot_mtx[:, i]
    b_col = b_mtx[:, i]

    den = (
        tot_col.toarray().ravel() if issparse(tot_col) else np.asarray(tot_col).ravel()
    )
    num = b_col.toarray().ravel() if issparse(b_col) else np.asarray(b_col).ravel()

    out = np.full_like(den, np.nan, dtype=np.float32)
    return np.divide(num, den, out=out, where=(den > 0))


def pseudobulk_by_groups(mat, group_idx, n_groups):
    """Sum the columns of ``mat`` within each group.

    Parameters
    ----------
    mat : sparse or ndarray
        Feature x cell matrix.
    group_idx : np.ndarray
        Length-n_cells int array; ``group_idx[c]`` is the group of column c.
    n_groups : int
        Number of groups (output column count).

    Returns
    -------
    np.ndarray
        Dense (n_features, n_groups) array of summed counts.
    """
    n_cells = mat.shape[1]
    assert group_idx.shape[0] == n_cells, "group_idx length must match #cells"
    indicator = csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (np.arange(n_cells), group_idx)),
        shape=(n_cells, n_groups),
    )
    out = mat @ indicator
    return out.toarray() if issparse(out) else np.asarray(out)


def compute_af_by_groups(tot_mtx, b_mtx, group_idx, n_groups):
    """Per-group pseudobulk allele frequency matrix of shape (n_features, n_groups)."""
    tot_grp = pseudobulk_by_groups(tot_mtx, group_idx, n_groups)
    b_grp = pseudobulk_by_groups(b_mtx, group_idx, n_groups)
    out = np.full(tot_grp.shape, np.nan, dtype=np.float32)
    return np.divide(b_grp, tot_grp, out=out, where=(tot_grp > 0))


def compute_af_pseudobulk(tot_mtx, b_mtx):
    """Compute per-SNP allele frequency across all cells (pseudobulk sum).

    Parameters
    ----------
    tot_mtx : sparse or ndarray
        Total depth matrix (SNPs x cells).
    b_mtx : sparse or ndarray
        B-allele count matrix.

    Returns
    -------
    np.ndarray
        Pseudobulk allele frequency per SNP; ``NaN`` where total depth is zero.
    """
    if issparse(tot_mtx):
        den = np.asarray(tot_mtx.sum(axis=1)).ravel()
    else:
        den = tot_mtx.sum(axis=1)

    if issparse(b_mtx):
        num = np.asarray(b_mtx.sum(axis=1)).ravel()
    else:
        num = b_mtx.sum(axis=1)

    out = np.full(den.shape[0], np.nan, dtype=np.float32)
    return np.divide(num, den, out=out, where=(den > 0))


##################################################
def get_mask_by_region(snps: pd.DataFrame, regions: pd.DataFrame) -> np.ndarray:
    """
    Return a boolean mask (len == len(snps)) indicating whether each SNP (CHR, POS)
    overlaps any interval in a BED-like file, using 0-based half-open intervals
    [Start, End).

    Assumes SNP POS is 1-based.
    """
    n = len(snps)
    r_chr = regions["#CHR"] if "#CHR" in regions.columns else regions["Chromosome"]
    r_start = (
        regions["START"] if "START" in regions.columns else regions["Start"]
    ).to_numpy()
    r_end = (regions["END"] if "END" in regions.columns else regions["End"]).to_numpy()

    keep = np.zeros(n, dtype=bool)
    for chrom in snps["#CHR"].unique():
        sm = (snps["#CHR"] == chrom).to_numpy()
        rm = (r_chr == chrom).to_numpy()
        if not rm.any():
            continue
        positions = snps.loc[sm, "POS"].to_numpy().astype(np.int64) - 1
        reg_starts = r_start[rm]
        reg_ends = r_end[rm]
        sort_idx = np.argsort(reg_starts)
        reg_starts, reg_ends = reg_starts[sort_idx], reg_ends[sort_idx]
        right_bounds = np.searchsorted(reg_starts, positions, side="right")
        overlaps = np.array(
            [
                np.any(reg_ends[: right_bounds[i]] > positions[i])
                for i in range(len(positions))
            ]
        )
        keep[np.where(sm)[0]] = overlaps
    return keep


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


##################################################
def subset_baf(
    baf_df: pd.DataFrame, ch: str, start: int, end: int, is_last_block=False
):
    """Slice a BAF DataFrame to a chromosomal interval ``[start, end)``.

    For the last block, the interval is closed on the right: ``[start, end]``.

    Parameters
    ----------
    baf_df : pd.DataFrame
        DataFrame with ``#CHR`` and ``POS`` columns (or ``POS`` as index).
    ch : str or None
        Chromosome to filter on; if None, no chromosome filter is applied.
    start, end : int
        Genomic position boundaries.
    is_last_block : bool
        If True, use a closed right boundary.

    Returns
    -------
    pd.DataFrame
        Filtered subset.
    """
    if ch != None:
        baf_ch = baf_df[baf_df["#CHR"] == ch]
    else:
        baf_ch = baf_df
    if baf_ch.index.name == "POS":
        pos = baf_ch.index
    else:
        pos = baf_ch["POS"]
    if is_last_block:
        return baf_ch[(pos >= start) & (pos <= end)]
    else:
        return baf_ch[(pos >= start) & (pos < end)]


def assign_snp_bounderies(
    snps: pd.DataFrame, regions: pd.DataFrame, colname="region_id"
):
    """
    divide regions into [START, END) subregions, each subregion has one SNP.
    If a SNP is out-of-region, its START and END will be 0 and region_id will be "".
    region_id is taken from regions["region_id"] (4th-column seg_id, with
    per-row fallback to "CHR:START-END" handled by read_region_file).
    """
    snps["START"] = 0
    snps["END"] = 0

    snps[colname] = ""

    chroms = snps["#CHR"].unique().tolist()
    region_grps_ch = regions.groupby(by="#CHR", sort=False)
    for chrom in chroms:
        regions_ch = region_grps_ch.get_group(chrom)
        for _, region in regions_ch.iterrows():
            reg_start, reg_end = region["START"], region["END"]
            reg_snps = subset_baf(snps, chrom, reg_start, reg_end)
            if len(reg_snps) == 0:
                continue
            reg_snp_positions = reg_snps["POS0"].to_numpy()
            reg_snp_indices = reg_snps.index.to_numpy()

            snps.loc[reg_snp_indices, colname] = region["region_id"]

            if len(reg_snps) == 1:
                snps.loc[reg_snp_indices, "START"] = reg_start
                snps.loc[reg_snp_indices, "END"] = reg_end
            else:
                reg_bounderies = np.ceil(
                    np.vstack([reg_snp_positions[:-1], reg_snp_positions[1:]]).mean(
                        axis=0
                    )
                ).astype(np.uint32)
                reg_bounderies = np.concatenate(
                    [[reg_start], reg_bounderies, [reg_end]]
                )
                snps.loc[reg_snp_indices, "START"] = reg_bounderies[:-1]
                snps.loc[reg_snp_indices, "END"] = reg_bounderies[1:]

    snps["BLOCKSIZE"] = snps["END"] - snps["START"]
    return snps
