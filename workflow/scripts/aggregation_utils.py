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


def count_split_genes(snps, grp_cols, bb_col="bb_id", feature_col="feature_id"):
    """Count genes whose SNPs cross a bin boundary (sanity check for gene-aware binning).

    A gene is "split" if, within one ``grp_cols`` group (region_id / PS / phase_group),
    its SNPs land in more than one bin. Returns ``(n_split, n_genes)`` over genic SNPs,
    or ``None`` if no ``feature_col`` is present.
    """
    if feature_col not in snps.columns:
        return None
    g = snps[snps[feature_col].notna() & (snps[feature_col] != "intergenic")]
    if len(g) == 0:
        return 0, 0
    key = g.groupby(grp_cols + [feature_col], sort=False)[bb_col].nunique()
    return int((key > 1).sum()), int(len(key))


def gene_block_labels(n_items, ranges):
    """Gene-block id per ordered item (windows or SNPs) so bins never split a gene.

    Each gene occupies an inclusive index range ``(lo, hi)`` over the item ordering
    (e.g. the first..last window holding that gene's SNPs). Every boundary internal to
    a range is made non-cuttable; overlapping ranges (genes that share an item, i.e.
    adjacent/overlapping genes) merge transitively into one block; items in no range
    are singleton blocks (native window/SNP granularity). Items sharing a block id must
    stay in one bin, so a bin holds whole genes only — never a partial gene.

    Parameters
    ----------
    n_items : int
        Number of ordered items.
    ranges : iterable of (lo, hi)
        Inclusive index ranges, one per gene.

    Returns
    -------
    np.ndarray (int64), length n_items
        Run-length-contiguous block id; the boundary between items ``i-1`` and ``i`` is
        a unit boundary iff ``labels[i] != labels[i-1]``.
    """
    if n_items == 0:
        return np.zeros(0, dtype=np.int64)
    blocked = np.zeros(n_items - 1, dtype=bool)  # blocked[i] = boundary (i, i+1) non-cuttable
    for lo, hi in ranges:
        if hi > lo:
            blocked[lo:hi] = True
    labels = np.empty(n_items, dtype=np.int64)
    labels[0] = 0
    if n_items > 1:
        labels[1:] = np.cumsum(~blocked)
    return labels


@numba.njit
def _bin_windows_numba(
    win_reads,
    win_nsnps,
    min_snp_reads,
    min_snp_per_block,
    win_starts,
    win_ends,
    max_blocksize,
    unit_ids,
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
    unit_ids : (W,) int64
        Gene-unit id per window; a bin may only close at a unit boundary
        (``unit_ids[i] != unit_ids[i-1]``), so a gene is never split across bins.
        The cap (``max_blocksize``) also only forces a cut at a unit boundary, so a
        single gene larger than the cap stays whole. Pass ``np.arange(W)`` for no
        constraint.

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
        unit_boundary = unit_ids[i] != unit_ids[i - 1]
        if ((meets_reads and acc_n >= min_snp_per_block) or exceeds_size) and unit_boundary:
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
    gene_aware=False,
):
    """Window-based adaptive binning: merge consecutive windows until SNP thresholds are met.

    Parameters
    ----------
    windows : pd.DataFrame
        Windows with ``#CHR``, ``START``, ``END``, ``win_idx``, and grouping columns.
        When ``gene_aware``, must also carry a ``gene_block`` column (see
        ``gene_block_labels``) so bins never split a gene.
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
        if gene_aware:
            grp_units = grp_wins["gene_block"].to_numpy()
        else:
            grp_units = np.arange(len(grp_idxs), dtype=np.int64)

        local_bin_ids, n_bins = _bin_windows_numba(
            grp_reads,
            grp_nsnps,
            min_snp_reads,
            min_snp_per_block,
            grp_starts,
            grp_ends,
            max_blocksize,
            grp_units,
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


def locate_atac_fragment_file(ranger_dir):
    """Return the 10x ATAC fragment file inside a cellranger dir, or None."""
    for fname in ("atac_fragments.tsv.gz", "fragments.tsv.gz"):
        fpath = os.path.join(ranger_dir, fname)
        if os.path.exists(fpath):
            return fpath
    return None


def rna_h5ad_to_bb(h5ad_file, barcodes, bb_df, num_bbs, assay_type):
    """Aggregate per-cell RNA counts (h5ad from ``process_rna_anndata``) into bb bins.

    Each RNA feature (gene) is assigned to the bb bin it overlaps most (``feature_to_blocks``
    -> largest overlap, same mapping as copytyping's ``cnv_segmentation``); its per-cell counts
    are summed into that bin. Cells are reordered to ``barcodes`` so the columns match that
    assay's ``bb.*allele.npz`` matrices.

    Parameters
    ----------
    h5ad_file : str
        AnnData (cells x genes) with ``var`` carrying ``#CHR``, ``START``, ``END``.
    barcodes : sequence of str
        Cell barcodes (``"{raw}_{rep}"``) in matrix-column order (that assay's allele columns).
    bb_df : pd.DataFrame
        Bins with ``#CHR``, ``START``, ``END`` (0-based half-open) and ``bb_id``.
    num_bbs : int
        Number of bins (output rows).

    Returns
    -------
    scipy.sparse.csr_matrix, shape ``(num_bbs, n_cells)``, dtype int32.
    """
    adata = sc.read_h5ad(h5ad_file)
    barcodes = np.asarray(barcodes, dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    if len(missing):
        raise ValueError(
            f"{len(missing)} barcodes missing from {h5ad_file}, e.g. {missing[:5]}"
        )
    adata = adata[barcodes, :].copy()
    adata = feature_to_blocks(adata, bb_df, assay_type, block_idx="bb_id", drop_cols=False)
    x_count = matrix_segmentation(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    return x_count.astype(np.int32)


def atac_fragments_to_bb(
    frag_files, reps, barcodes_full, bb_df, num_bbs, chunksize=5_000_000
):
    """Count deduped ATAC fragments per bb bin per cell from 10x fragment files.

    Each row of a 10x ``atac_fragments.tsv.gz`` is one deduplicated fragment
    (``chrom, start, end, barcode, readSupport``); the readSupport column is IGNORED.
    Every fragment is counted once, assigned to the bb bin containing its midpoint, so
    the column sums equal the number of in-bin fragments per cell.

    Parameters
    ----------
    frag_files, reps : parallel lists
        ``frag_files[i]`` is the fragment file for replicate ``reps[i]``.
    barcodes_full : pd.DataFrame
        Columns ``REP_ID``, ``BARCODE`` (``BARCODE`` = ``"{raw}_{rep}"``) giving the cell
        column order (identical to that assay's ``bb.*allele.npz`` columns).
    bb_df : pd.DataFrame
        Bins with ``#CHR``, ``START``, ``END`` (0-based half-open) and ``bb_id``.
    num_bbs : int
        Number of bins (output rows).

    Returns
    -------
    scipy.sparse.csr_matrix, shape ``(num_bbs, n_cells)``, dtype int32.
    """
    n_cells = len(barcodes_full)
    bc_rep = barcodes_full["REP_ID"].to_numpy().astype(str)
    bc_full = barcodes_full["BARCODE"].to_numpy().astype(str)
    # global column index keyed by (rep, raw_barcode); strip the "_{rep}" suffix
    col_of = {}
    for i in range(n_cells):
        rep, raw = bc_rep[i], bc_full[i]
        sfx = "_" + rep
        if raw.endswith(sfx):
            raw = raw[: -len(sfx)]
        col_of[(rep, raw)] = i

    rows_all, cols_all = [], []
    for frag_file, rep in zip(frag_files, reps):
        rep_map = {raw: c for (r, raw), c in col_of.items() if r == rep}
        if not rep_map or frag_file is None:
            continue
        n_frag = 0
        for chunk in pd.read_csv(
            frag_file, sep="\t", comment="#", header=None, usecols=[0, 1, 2, 3],
            names=["#CHR", "start", "end", "BC"],
            dtype={0: str, 1: np.int64, 2: np.int64, 3: str}, chunksize=chunksize,
        ):
            col_vals = chunk["BC"].map(rep_map).to_numpy()
            m = ~pd.isna(col_vals)
            if not m.any():
                continue
            sub = chunk.loc[m]
            mid = (sub["start"].to_numpy() + sub["end"].to_numpy()) // 2
            frag = pd.DataFrame({"#CHR": sub["#CHR"].to_numpy(), "POS0": mid})
            frag = assign_pos_to_range(frag, bb_df, ref_id="bb_id", pos_col="POS0")
            keep = frag["bb_id"].notna().to_numpy()
            if not keep.any():
                continue
            rows_all.append(frag.loc[keep, "bb_id"].to_numpy().astype(np.int64))
            cols_all.append(col_vals[m][keep].astype(np.int64))
            n_frag += int(keep.sum())
        logging.info(f"  ATAC {rep}: {n_frag} in-bin fragments counted")

    if rows_all:
        rows = np.concatenate(rows_all)
        cols = np.concatenate(cols_all)
    else:
        rows = np.zeros(0, dtype=np.int64)
        cols = np.zeros(0, dtype=np.int64)
    data = np.ones(len(rows), dtype=np.int32)
    return csr_matrix((data, (rows, cols)), shape=(num_bbs, n_cells), dtype=np.int32)
