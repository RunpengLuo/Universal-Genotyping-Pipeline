import os
import logging

import numpy as np
import pandas as pd
import numba

from scipy.sparse import csr_matrix, hstack, issparse
from scipy.stats import beta as beta_dist

import scanpy as sc

from io_utils import *
from combine_counts_utils import *
from count_reads_utils import *


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


@numba.njit
def _bin_snps_numba(read_counts, min_snp_reads, min_snp_per_block):
    """Greedy adaptive binning over a contiguous read-count array.

    Parameters
    ----------
    read_counts : (N, M) contiguous float64
        Per-SNP read counts for tumor samples.
    min_snp_reads : int
        Minimum total reads per sample for a bin to be complete.
    min_snp_per_block : int
        Minimum number of SNPs per bin.

    Returns
    -------
    bin_ids : (N,) int64
        Relative bin ID for each SNP within this group.
    n_bins : int
        Number of bins created (before last-block adjustment).
    """
    N, M = read_counts.shape
    bin_ids = np.zeros(N, dtype=np.int64)
    if N == 0:
        return bin_ids, 0

    bin_id = 0
    prev_start = 0
    acc = read_counts[0].copy()
    acc_n = 1

    for i in range(1, N):
        meets_reads = True
        if min_snp_reads > 0:
            for j in range(M):
                if acc[j] < min_snp_reads:
                    meets_reads = False
                    break
        if meets_reads and acc_n >= min_snp_per_block:
            bin_ids[prev_start:i] = bin_id
            bin_id += 1
            prev_start = i
            acc = read_counts[i].copy()
            acc_n = 1
        else:
            for j in range(M):
                acc[j] += read_counts[i, j]
            acc_n += 1

    # last block: keep as own bin if it meets all criteria and isn't the only block
    last_meets_reads = True
    if min_snp_reads > 0:
        for j in range(M):
            if acc[j] < min_snp_reads:
                last_meets_reads = False
                break
    if last_meets_reads and acc_n >= min_snp_per_block and prev_start > 0:
        bin_ids[prev_start:] = bin_id
        bin_id += 1
    else:
        merge_id = bin_id - 1 if bin_id > 0 else 0
        bin_ids[prev_start:] = merge_id

    return bin_ids, bin_id


@numba.njit
def _bin_windows_numba(
    win_reads,
    win_nsnps,
    min_snp_reads,
    min_snp_per_block,
    win_starts,
    win_ends,
    max_blocksize,
):
    """Greedy adaptive binning over consecutive windows.

    Parameters
    ----------
    win_reads : (W, M) contiguous float64
        Per-window total tumor reads (summed from SNPs in each window).
    win_nsnps : (W,) int64
        Number of SNPs per window.
    min_snp_reads : int
        Minimum total reads per sample for a bin to be complete.
    min_snp_per_block : int
        Minimum number of SNPs per bin.
    win_starts : (W,) int64
        START coordinate per window.
    win_ends : (W,) int64
        END coordinate per window.
    max_blocksize : int
        Maximum genomic span (END - START) for a bin. When exceeded, force
        a bin boundary. Set to 0 to disable.

    Returns
    -------
    bin_ids : (W,) int64
        Relative bin ID for each window within this group.
    n_bins : int
        Number of bins created (before last-block adjustment).
    """
    W, M = win_reads.shape
    bin_ids = np.zeros(W, dtype=np.int64)
    if W == 0:
        return bin_ids, 0

    bin_id = 0
    prev_start = 0
    acc = win_reads[0].copy()
    acc_n = win_nsnps[0]

    for i in range(1, W):
        meets_reads = True
        if min_snp_reads > 0:
            for j in range(M):
                if acc[j] < min_snp_reads:
                    meets_reads = False
                    break
        span = win_ends[i - 1] - win_starts[prev_start]
        exceeds_size = max_blocksize > 0 and span >= max_blocksize
        if (meets_reads and acc_n >= min_snp_per_block) or exceeds_size:
            bin_ids[prev_start:i] = bin_id
            bin_id += 1
            prev_start = i
            acc = win_reads[i].copy()
            acc_n = win_nsnps[i]
        else:
            for j in range(M):
                acc[j] += win_reads[i, j]
            acc_n += win_nsnps[i]

    # last block: keep as own bin if it meets all criteria and isn't the only block
    last_meets_reads = True
    if min_snp_reads > 0:
        for j in range(M):
            if acc[j] < min_snp_reads:
                last_meets_reads = False
                break
    last_span = win_ends[W - 1] - win_starts[prev_start]
    last_exceeds_size = max_blocksize > 0 and last_span >= max_blocksize
    if (last_meets_reads and acc_n >= min_snp_per_block and prev_start > 0) or (
        last_exceeds_size and prev_start > 0
    ):
        bin_ids[prev_start:] = bin_id
        bin_id += 1
    else:
        merge_id = bin_id - 1 if bin_id > 0 else 0
        bin_ids[prev_start:] = merge_id

    return bin_ids, bin_id


def adaptive_segmentation(
    windows: pd.DataFrame,
    snps: pd.DataFrame,
    tot_mtx: np.ndarray,
    min_snp_reads: int,
    min_snp_per_block: int,
    grp_cols: list,
    tumor_sidx=0,
    max_blocksize=0,
):
    """Window-based adaptive binning: merge consecutive windows until SNP thresholds are met.

    Parameters
    ----------
    windows : pd.DataFrame
        Windows with ``#CHR``, ``START``, ``END``, ``win_idx``, and grouping columns.
    snps : pd.DataFrame
        SNP DataFrame with ``POS0`` and ``#CHR`` columns.
    tot_mtx : (N, M) ndarray
        Per-SNP total read counts (N SNPs, M samples).
    min_snp_reads : int
        Minimum total tumor reads per sample for a bin.
    min_snp_per_block : int
        Minimum number of SNPs per bin.
    grp_cols : list of str
        Columns to group windows by (e.g. ``["region_id"]``).
    tumor_sidx : int
        Index of first tumor sample column.

    Returns
    -------
    bbs : pd.DataFrame
        Bin definitions with ``#CHR``, ``START``, ``END``, ``#SNPS``, ``BLOCKSIZE``,
        ``bb_id``, and grouping columns.
    snps : pd.DataFrame
        Input SNPs with ``bb_id`` and ``win_idx`` columns added. SNPs not falling
        in any window are dropped.
    """
    from scipy.sparse import issparse

    logging.info(
        f"adaptive_segmentation: min_snp_reads={min_snp_reads}, "
        f"min_snp_per_block={min_snp_per_block}, "
        f"max_blocksize={max_blocksize}"
    )

    # 1. Assign SNPs to windows
    snps["_orig_idx"] = np.arange(len(snps))
    snps = assign_pos_to_range(snps, windows, ref_id="win_idx", pos_col="POS0")
    outside_mask = snps["win_idx"].isna()
    n_outside = outside_mask.sum()
    logging.info(
        f"SNPs outside any window: {n_outside}/{len(snps)} ({n_outside / max(len(snps), 1):.3%})"
    )
    if n_outside > 0:
        off_idx = snps.loc[outside_mask, "_orig_idx"].to_numpy()
        off_depth = tot_mtx[off_idx].sum(axis=1)
        logging.info(
            f"off-target SNP depth: "
            f"min={off_depth.min()}, max={off_depth.max()}, "
            f"mean={off_depth.mean():.1f}, median={np.median(off_depth):.1f}"
        )
    snps = snps.dropna(subset=["win_idx"]).reset_index(drop=True)
    snps["win_idx"] = snps["win_idx"].astype(np.int64)

    # 2. Compute per-window stats
    W = len(windows)
    M_tumor = tot_mtx.shape[1] - tumor_sidx
    win_nsnps = np.zeros(W, dtype=np.int64)
    win_reads = np.zeros((W, M_tumor), dtype=np.float64)

    snp_win_idx = snps["win_idx"].to_numpy()
    snp_orig_idx = snps["_orig_idx"].to_numpy()

    if issparse(tot_mtx):
        tot_tumor = tot_mtx[:, tumor_sidx:].toarray().astype(np.float64)
    else:
        tot_tumor = tot_mtx[:, tumor_sidx:].astype(np.float64)

    for i in range(len(snps)):
        w = snp_win_idx[i]
        win_nsnps[w] += 1
        win_reads[w] += tot_tumor[snp_orig_idx[i]]

    # 3. Group windows and run numba kernel
    bin_id = 0
    windows["bin_id"] = 0
    all_win_starts = windows["START"].to_numpy(dtype=np.int64)
    all_win_ends = windows["END"].to_numpy(dtype=np.int64)
    win_grps = windows.groupby(by=grp_cols, sort=False)
    logging.info(f"#window groups={len(win_grps)}, grouper: {grp_cols}")

    for _, grp_wins in win_grps:
        grp_idxs = grp_wins.index.to_numpy()
        grp_reads = np.ascontiguousarray(win_reads[grp_idxs])
        grp_nsnps = np.ascontiguousarray(win_nsnps[grp_idxs])
        grp_starts = np.ascontiguousarray(all_win_starts[grp_idxs])
        grp_ends = np.ascontiguousarray(all_win_ends[grp_idxs])

        local_bin_ids, n_bins = _bin_windows_numba(
            grp_reads,
            grp_nsnps,
            min_snp_reads,
            min_snp_per_block,
            grp_starts,
            grp_ends,
            max_blocksize,
        )

        local_bin_ids += bin_id
        windows.loc[grp_idxs, "bin_id"] = local_bin_ids
        bin_id += max(n_bins, 1)

    # 4. Propagate bin_id to SNPs
    win_bin_map = windows["bin_id"].to_numpy()
    snps["bb_id"] = win_bin_map[snps["win_idx"].to_numpy()]

    # 5. Build bbs DataFrame from windows
    pos_dict = {
        "#CHR": ("#CHR", "first"),
        "START": ("START", "min"),
        "END": ("END", "max"),
    }
    for grp_col in grp_cols:
        pos_dict[grp_col] = (grp_col, "first")

    win_grps_by_bin = windows.groupby("bin_id", sort=True)
    bbs = win_grps_by_bin.agg(**pos_dict)

    # SNP counts per bin
    snp_counts = snps.groupby("bb_id").size()
    bbs["#SNPS"] = bbs.index.map(snp_counts).fillna(0).astype(int)
    bbs["BLOCKSIZE"] = bbs["END"] - bbs["START"]
    bbs["bb_id"] = bbs.index

    # PS column if present
    if "PS" in snps.columns:
        ps = snps.groupby("bb_id")["PS"].first()
        bbs["PS"] = bbs.index.map(ps)

    num_bbs = len(bbs)
    bin_sizes = bbs["#SNPS"].to_numpy()
    block_sizes = bbs["BLOCKSIZE"].to_numpy()
    logging.info("adaptive_segmentation summary")
    logging.info(f"#SNPs={len(snps)}, #windows={W}, #bins={num_bbs}")
    logging.info(
        "snps per bin: min=%.0f  median=%.0f  max=%.0f",
        float(bin_sizes.min()) if num_bbs > 0 else 0,
        float(np.median(bin_sizes)) if num_bbs > 0 else 0,
        float(bin_sizes.max()) if num_bbs > 0 else 0,
    )
    logging.info(
        "blocksize per bin: min=%.0f  median=%.0f  max=%.0f",
        float(block_sizes.min()) if num_bbs > 0 else 0,
        float(np.median(block_sizes)) if num_bbs > 0 else 0,
        float(block_sizes.max()) if num_bbs > 0 else 0,
    )

    return bbs, snps


def _assign_chrom_overlapping(qry, qry_mask, ref_chrom, ref_id, pos_col):
    """Assign positions to overlapping intervals on one chromosome using numpy."""
    positions = qry.loc[qry_mask, pos_col].to_numpy()
    qry_indices = qry.index[qry_mask]
    starts = ref_chrom["START"].to_numpy()
    ends = ref_chrom["END"].to_numpy()
    ids = ref_chrom[ref_id].to_numpy()

    sort_idx = np.argsort(starts)
    starts = starts[sort_idx]
    ends = ends[sort_idx]
    ids = ids[sort_idx]

    right_bounds = np.searchsorted(starts, positions, side="right")
    for i in range(len(positions)):
        pos = positions[i]
        cands = slice(0, right_bounds[i])
        mask = ends[cands] > pos
        if mask.any():
            qry.loc[qry_indices[i], ref_id] = ids[cands][mask][0]


def assign_pos_to_range(
    qry: pd.DataFrame,
    ref: pd.DataFrame,
    ref_id="region_id",
    pos_col="POS0",
    nodup=True,
):
    """Assign each query position to the reference interval it falls within.

    Uses ``np.searchsorted`` for non-overlapping intervals (fast path) and
    falls back to a loop for chromosomes with overlapping intervals.

    Parameters
    ----------
    qry : pd.DataFrame
        Query DataFrame with ``#CHR`` and *pos_col* columns.
    ref : pd.DataFrame
        Reference intervals with ``#CHR``, ``START``, ``END``, and *ref_id*.
    ref_id : str
        Column name for the reference interval identifier.
    pos_col : str
        Column in *qry* containing 0-based positions.
    nodup : bool
        If True, each query is assigned to at most one reference interval;
        if False, returns all overlapping (query, ref) pairs.

    Returns
    -------
    pd.DataFrame
        *qry* with *ref_id* column added (if ``nodup=True``), or a hits
        DataFrame (if ``nodup=False``).
    """
    if not nodup:
        rows = []
        for chrom in ref["#CHR"].unique():
            qm = (qry["#CHR"] == chrom).to_numpy()
            if not qm.any():
                continue
            positions = qry.loc[qm, pos_col].to_numpy()
            q_indices = qry.index[qm].to_numpy()
            starts = ref.loc[ref["#CHR"] == chrom, "START"].to_numpy()
            ends = ref.loc[ref["#CHR"] == chrom, "END"].to_numpy()
            ids = ref.loc[ref["#CHR"] == chrom, ref_id].to_numpy()
            sort_idx = np.argsort(starts)
            starts, ends, ids = starts[sort_idx], ends[sort_idx], ids[sort_idx]
            right_bounds = np.searchsorted(starts, positions, side="right")
            for i in range(len(positions)):
                pos = positions[i]
                cands = slice(0, right_bounds[i])
                mask = ends[cands] > pos
                for rid in ids[cands][mask]:
                    rows.append((chrom, pos, q_indices[i], rid))
        hits = pd.DataFrame(rows, columns=["#CHR", "POS0", "qry_index", ref_id])
        hits["POS"] = hits["POS0"] + 1
        return hits

    # nodup=True: assign each query position to at most one interval
    qry[ref_id] = pd.NA
    for chrom in ref["#CHR"].unique():
        qry_mask = (qry["#CHR"] == chrom).to_numpy()
        if not qry_mask.any():
            continue

        ref_chrom = ref.loc[ref["#CHR"] == chrom].sort_values("START")
        starts = ref_chrom["START"].to_numpy()
        ends = ref_chrom["END"].to_numpy()
        ids = ref_chrom[ref_id].to_numpy()
        positions = qry.loc[qry_mask, pos_col].to_numpy()

        has_overlap = len(starts) > 1 and np.any(starts[1:] < ends[:-1])

        if not has_overlap:
            # Fast path: searchsorted for non-overlapping intervals
            # 0-based half-open: START <= pos < END
            idx = np.searchsorted(starts, positions, side="right") - 1
            safe_idx = idx.clip(min=0)
            valid = (idx >= 0) & (positions < ends[safe_idx])
            qry_indices = qry.index[qry_mask]
            qry.loc[qry_indices[valid], ref_id] = ids[idx[valid]]
        else:
            _assign_chrom_overlapping(qry, qry_mask, ref_chrom, ref_id, pos_col)

    return qry


def snp_to_region(
    snp_df: pd.DataFrame, region_df: pd.DataFrame, assay_type: str, region_id="BIN_ID"
):
    """
    region_df must be 0-indexed non-overlapping intervals [s, t) in standard BED format.
    """
    logging.info(f"#{assay_type}-SNP (raw)={len(snp_df)}")
    snp_df = assign_pos_to_range(snp_df, region_df, ref_id=region_id, pos_col="POS0")
    isna_snp_df = snp_df[region_id].isna()
    logging.info(
        f"#{assay_type}: #SNPS outside any region={np.sum(isna_snp_df) / len(snp_df):.3%}"
    )
    snp_df.dropna(subset=region_id, inplace=True)
    snp_df[region_id] = snp_df[region_id].astype(region_df[region_id].dtype)
    logging.info(f"#{assay_type}-SNP (remain)={len(snp_df)}")

    counts = snp_df[region_id].value_counts()
    region_df[f"#SNPS"] = region_df[region_id].map(counts).fillna(0).astype(int)
    return snp_df


def adaptive_binning(
    snps: pd.DataFrame,
    min_snp_reads: int,
    min_snp_per_block: int,
    tot_mtx: np.ndarray,
    grp_cols: list,
    colname="",
    tumor_sidx=0,
):
    """
    Adaptive binning over pre-defined chromosome regions.
    Each bin satisfies MSR (min total reads) and MSPB (min SNPs per block).
    The last block in each group is kept as its own bin if it meets all criteria,
    otherwise it is merged into the previous bin.
    In-place adds <colname> column to snps to indicate bin ids.
    """
    bin_id = 0
    snps[colname] = 0
    snp_grps = snps.groupby(by=grp_cols, sort=False)
    logging.info(
        f"adaptive binning, colname={colname}, min_snp_reads={min_snp_reads}, min_snp_per_block={min_snp_per_block}"
    )
    logging.info(f"#SNP groups={len(snp_grps)}, grouper: {grp_cols}")
    for _, grp_snps in snp_grps:
        grp_idxs = grp_snps.index.to_numpy()

        grp_mtx = tot_mtx[grp_idxs]
        if issparse(grp_mtx):
            grp_reads = np.ascontiguousarray(
                grp_mtx[:, tumor_sidx:].toarray(), dtype=np.float64
            )
        else:
            grp_reads = np.ascontiguousarray(
                grp_mtx[:, tumor_sidx:], dtype=np.float64
            ).copy()

        local_bin_ids, n_bins = _bin_snps_numba(
            grp_reads, min_snp_reads, min_snp_per_block
        )

        local_bin_ids += bin_id
        snps.loc[grp_idxs, colname] = local_bin_ids
        bin_id += max(n_bins, 1)

    pos_dict = {
        "#CHR": ("#CHR", "first"),
        "START": ("START", "min"),
        "END": ("END", "max"),
        "START0": ("POS0", "min"),
        "END0": ("POS", "max"),
    }
    for grp_col in grp_cols:
        pos_dict[grp_col] = (grp_col, "first")

    snp_grps = snps.groupby(by=colname, sort=False, as_index=True)
    snp_bins = snp_grps.agg(**pos_dict)
    for feat_col in ("feature_id", "feature_type"):
        if feat_col in snps.columns:
            snp_bins[feat_col] = snp_grps[feat_col].agg(
                lambda x: ";".join(dict.fromkeys(str(v) for v in x if pd.notna(v)))
            )
    snp_bins.loc[:, "#SNPS"] = snp_grps.size()  # align by bin_id index, not position
    snp_bins.loc[:, "BLOCKSIZE"] = snp_bins["END"] - snp_bins["START"]
    snp_bins = snp_bins.reset_index()

    bin_sizes = snp_bins["#SNPS"].to_numpy()
    block_sizes = snp_bins["BLOCKSIZE"].to_numpy()
    logging.info("adaptive binning summary")
    logging.info(f"#SNPs={len(snps)}")
    logging.info(f"#bins={len(snp_bins)}")
    logging.info(
        "snps per bin: min=%.0f  median=%.0f  max=%.0f",
        float(bin_sizes.min()),
        float(np.median(bin_sizes)),
        float(bin_sizes.max()),
    )
    logging.info(
        "blocksize per bin: min=%.0f  median=%.0f  max=%.0f",
        float(block_sizes.min()),
        float(np.median(block_sizes)),
        float(block_sizes.max()),
    )

    return snp_bins


def matrix_segmentation(X, bin_ids, K):
    """
    N: #features
    M: #samples
    X: (N, M) sparse or dense   [snp-by-sample]
    bin_ids: (N,) ints in [0..K-1]  (assign each SNP to a bin, bin ids are non-decreasing)
    return:
      - sparse in  -> (K, M) csr_matrix
      - dense in   -> (K, M) ndarray
    """
    X = X.tocsr() if issparse(X) else np.asarray(X)

    bin_ids = np.asarray(bin_ids, dtype=np.int64)
    N, M = X.shape
    if bin_ids.shape[0] != N:
        raise ValueError(f"bin_ids length {bin_ids.shape[0]} != N {N}")
    if N and (bin_ids.min() < 0 or bin_ids.max() >= K):
        raise ValueError("bin_ids out of range")

    # (K, N) one-hot: row=bin, col=snp
    B = csr_matrix(
        (np.ones(N, dtype=np.int8), (bin_ids, np.arange(N, dtype=np.int64))),
        shape=(K, N),
    )

    X_bin = B @ X  # (K, M)
    return X_bin


def assign_largest_overlap(
    qry: pd.DataFrame, ref: pd.DataFrame, qry_id: str, ref_id: str
) -> pd.DataFrame:
    """For each row in qry, assign the ID of the overlapping ref interval
    with the largest overlap length.

    Both *qry* and *ref* must have ``#CHR``, ``START``, ``END`` columns
    (0-based half-open).
    """
    qry = qry.copy()
    qry[ref_id] = pd.NA

    for chrom in ref["#CHR"].unique():
        qm = qry["#CHR"] == chrom
        rm = ref["#CHR"] == chrom
        if not qm.any():
            continue

        q_starts = qry.loc[qm, "START"].to_numpy()
        q_ends = qry.loc[qm, "END"].to_numpy()
        r_starts = ref.loc[rm, "START"].to_numpy()
        r_ends = ref.loc[rm, "END"].to_numpy()
        r_ids = ref.loc[rm, ref_id].to_numpy()

        sort_idx = np.argsort(r_starts)
        r_starts = r_starts[sort_idx]
        r_ends = r_ends[sort_idx]
        r_ids = r_ids[sort_idx]

        right_bounds = np.searchsorted(r_starts, q_ends, side="left")

        best_ids = np.empty(len(q_starts), dtype=object)
        best_ids[:] = pd.NA
        for i in range(len(q_starts)):
            qs, qe = q_starts[i], q_ends[i]
            cands = slice(0, right_bounds[i])
            mask = r_ends[cands] > qs
            if not mask.any():
                continue
            c_starts = r_starts[cands][mask]
            c_ends = r_ends[cands][mask]
            c_ids = r_ids[cands][mask]
            overlap = np.minimum(qe, c_ends) - np.maximum(qs, c_starts)
            best_ids[i] = c_ids[np.argmax(overlap)]

        qry.loc[qm, ref_id] = best_ids

    return qry


def feature_to_blocks(
    adata: sc.AnnData,
    blocks: pd.DataFrame,
    assay_type: str,
    feature_idx="feature_idx",
    block_idx="region_id",
    drop_cols=True,
):
    """
    filter features not in blocks, likely masked regions include centromeres
    """
    logging.info(f"assign {assay_type} features to blocks, {feature_idx}-{block_idx}")
    if adata.is_view:
        adata = adata.copy()
    adata.var[feature_idx] = np.arange(len(adata.var))

    feature_df = adata.var.reset_index(drop=True)
    logging.info(f"#{assay_type}-features (raw)={len(feature_df)}")

    feature_df = assign_largest_overlap(feature_df, blocks, feature_idx, block_idx)
    isna_features = feature_df[block_idx].isna()
    logging.info(
        f"#{assay_type} feature outside any blocks={np.sum(isna_features) / len(feature_df):.3%}"
    )
    feature_df.dropna(subset=block_idx, inplace=True)
    feature_df[block_idx] = feature_df[block_idx].astype(blocks[block_idx].dtype)
    logging.info(f"#{assay_type} feature (remain)={len(feature_df)}")

    ##################################################
    adata.var = (
        adata.var.reset_index(drop=False)
        .merge(
            right=feature_df[[feature_idx, block_idx]],
            on=feature_idx,
            how="left",
        )
        .set_index("index")
    )
    adata = adata[:, adata.var[block_idx].notna()].copy()
    if drop_cols:
        adata.var.drop(columns=[feature_idx, block_idx], inplace=True)
    else:
        adata.var[block_idx] = adata.var[block_idx].astype(feature_df[block_idx].dtype)
    return adata
