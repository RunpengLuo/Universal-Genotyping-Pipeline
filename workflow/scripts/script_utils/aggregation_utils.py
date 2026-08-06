import logging

import numpy as np
import pandas as pd
import numba

from scipy.sparse import issparse

from interval_utils import (
    assign_all_features,
    assign_pos_to_range,
    overlaps_any_range,
)
from io_utils import read_BED, read_gtf
from matrix_utils import group_sum


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
    blocked = np.zeros(
        n_items - 1, dtype=bool
    )  # blocked[i] = boundary (i, i+1) non-cuttable
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
    min_snp_reads_vec,
    min_snp_per_bin,
    win_starts,
    win_ends,
    max_blocksize,
    unit_ids,
):
    """Greedy adaptive binning over consecutive windows.

    A bin can close only when its accumulated reads meet the per-column threshold
    (``acc[j] >= min_snp_reads_vec[j]`` for every sample j). ``max_blocksize`` is a
    SECONDARY guard: it may force a cut only once that read threshold is already
    satisfied, so no bin is ever emitted below threshold and every bin is non-empty
    in every column (e.g. every WES column carries reads).

    Parameters
    ----------
    win_reads : (W, M) contiguous float64
        Per-window total tumor reads (summed from SNPs in each window).
    win_nsnps : (W,) int64
        Number of SNPs per window.
    min_snp_reads_vec : (M,) float64
        Per-column minimum SNP reads per bin.
    min_snp_per_bin : int
        Minimum number of SNPs per bin.
    win_starts : (W,) int64
        START coordinate per window.
    win_ends : (W,) int64
        END coordinate per window.
    max_blocksize : int
        Span cap (bp); once the read threshold is met, a bin over this span is cut
        even if it holds fewer than ``min_snp_per_bin`` SNPs. Never cuts below the
        read threshold. Set to 0 to disable.
    unit_ids : (W,) int64
        Gene-unit id per window; a bin may only close at a unit boundary
        (``unit_ids[i] != unit_ids[i-1]``), so a gene is never split across bins.
        Pass ``np.arange(W)`` for no constraint.

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
        for j in range(M):
            if acc[j] < min_snp_reads_vec[j]:
                meets_reads = False
                break
        span = win_ends[i - 1] - win_starts[prev_start]
        exceeds_size = max_blocksize > 0 and span >= max_blocksize
        unit_boundary = unit_ids[i] != unit_ids[i - 1]
        if meets_reads and (acc_n >= min_snp_per_bin or exceeds_size) and unit_boundary:
            bin_ids[prev_start:i] = bin_id
            bin_id += 1
            prev_start = i
            acc = win_reads[i].copy()
            acc_n = win_nsnps[i]
        else:
            for j in range(M):
                acc[j] += win_reads[i, j]
            acc_n += win_nsnps[i]

    # last block: keep as own bin if it meets the read threshold and isn't the only block
    last_meets_reads = True
    for j in range(M):
        if acc[j] < min_snp_reads_vec[j]:
            last_meets_reads = False
            break
    last_span = win_ends[W - 1] - win_starts[prev_start]
    last_exceeds_size = max_blocksize > 0 and last_span >= max_blocksize
    if (
        last_meets_reads
        and (acc_n >= min_snp_per_bin or last_exceeds_size)
        and prev_start > 0
    ):
        bin_ids[prev_start:] = bin_id
        bin_id += 1
    else:
        merge_id = bin_id - 1 if bin_id > 0 else 0
        bin_ids[prev_start:] = merge_id

    return bin_ids, bin_id


def assign_and_keep(snps, ref, ref_id, pos_col="POS0", label=""):
    """Assign SNPs to reference intervals, log the misses, and drop them.

    The shared body of ``snps_to_windows`` and ``snp_to_region``: assign, report the
    fraction outside every interval, drop those rows, and cast the id to the
    reference's own dtype.

    Args:
        snps: SNP frame with ``#CHR`` and *pos_col*.
        ref: Reference intervals carrying *ref_id*.
        ref_id: Identifier column to assign.
        pos_col: 0-based position column of *snps*.
        label: Prefix for the log line.

    Returns:
        ``(kept, outside_mask)``: the assigned SNPs reindexed from 0, and the
        boolean mask of dropped rows over the INPUT frame.
    """
    snps = assign_pos_to_range(snps, ref, ref_id=ref_id, pos_col=pos_col)
    outside = snps[ref_id].isna()
    n_out = int(outside.sum())
    logging.info(
        f"{label}SNPs outside any {ref_id}: {n_out}/{len(snps)} "
        f"({n_out / max(len(snps), 1):.3%})"
    )
    kept = snps.loc[~outside].reset_index(drop=True)
    kept[ref_id] = kept[ref_id].astype(ref[ref_id].dtype)
    return kept, outside.to_numpy()


def snps_to_windows(snps, windows, tot_mtx):
    """Assign SNPs to windows once, dropping those outside every window.

    Hoisted out of ``adaptive_segmentation`` so a caller sweeping a binning grid
    pays for the assignment once instead of per grid point. The returned frame
    carries ``_orig_idx`` (row in the input frame, indexing ``tot_mtx``) and an
    int64 ``win_idx``; ``adaptive_segmentation`` skips step 1 when it sees one.

    Args:
        snps: SNP frame with ``#CHR`` and ``POS0``.
        windows: Window frame with ``#CHR``, ``START``, ``END``, ``win_idx``.
        tot_mtx: Per-SNP total counts, for the off-window depth log only.

    Returns:
        The SNPs inside a window, reindexed from 0.
    """
    snps["_orig_idx"] = np.arange(len(snps))
    orig_idx = snps["_orig_idx"].to_numpy()
    kept, outside = assign_and_keep(snps, windows, "win_idx")
    if outside.any():
        off_depth = tot_mtx[orig_idx[outside]].sum(axis=1)
        logging.info(
            f"off-target SNP depth: "
            f"min={off_depth.min()}, max={off_depth.max()}, "
            f"mean={off_depth.mean():.1f}, median={np.median(off_depth):.1f}"
        )
    kept["win_idx"] = kept["win_idx"].astype(np.int64)
    return kept


def adaptive_segmentation(
    windows: pd.DataFrame,
    snps: pd.DataFrame,
    tot_mtx: np.ndarray,
    min_snp_reads,
    min_snp_per_bin: int,
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
    min_snp_reads : int or array-like
        Minimum total tumor reads for a bin, per tumor column. A scalar is
        broadcast to every tumor column; an array of length ``M - tumor_sidx``
        sets a per-column threshold.
    min_snp_per_bin : int
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

    M_tumor = tot_mtx.shape[1] - tumor_sidx
    min_snp_reads_vec = np.ascontiguousarray(
        np.broadcast_to(np.asarray(min_snp_reads, dtype=np.float64), (M_tumor,))
    )
    logging.info(
        f"adaptive_segmentation: min_snp_reads(per-col)={min_snp_reads_vec.tolist()}, "
        f"min_snp_per_bin={min_snp_per_bin}, "
        f"max_blocksize={max_blocksize}"
    )

    # 1. Assign SNPs to windows, unless the caller already did (see snps_to_windows)
    if "win_idx" not in snps.columns:
        snps = snps_to_windows(snps, windows, tot_mtx)

    # 2. Compute per-window stats
    W = len(windows)
    snp_win_idx = snps["win_idx"].to_numpy()
    snp_orig_idx = snps["_orig_idx"].to_numpy()

    if issparse(tot_mtx):
        tot_tumor = tot_mtx[:, tumor_sidx:].toarray().astype(np.float64)
    else:
        tot_tumor = tot_mtx[:, tumor_sidx:].astype(np.float64)

    win_nsnps = np.bincount(snp_win_idx, minlength=W).astype(np.int64)
    win_reads = np.asarray(
        group_sum(tot_tumor[snp_orig_idx], snp_win_idx, W, axis=0), dtype=np.float64
    )

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
            min_snp_reads_vec,
            min_snp_per_bin,
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


def merge_feature_ids(strings, sep=";", default="intergenic"):
    """Collapse an iterable of ``sep``-joined feature_id strings into one deduped union.

    Drops ``default`` tokens unless nothing else remains; preserves first-seen order.
    """
    seen = dict()
    for s in strings:
        if not isinstance(s, str):
            continue
        for tok in s.split(sep):
            if tok and tok != default:
                seen[tok] = None
    return sep.join(seen) if seen else default


def annotate_feature_type(snps, gtf_file):
    """Annotate SNPs with feature_id (overlapping genes) and feature_type.

    ``feature_id`` is a ``;``-joined list of every GTF gene the SNP overlaps
    (``intergenic`` when none), so gene membership is read back off it rather than
    assigned a second time. ``feature_type`` is exon > intron > intergenic.
    """
    gtf = read_gtf(gtf_file, ("gene", "exon"))
    snps["feature_id"] = assign_all_features(snps, gtf["gene"], id_col="gene_id")
    in_gene = snps["feature_id"] != "intergenic"
    in_exon = pd.Series(overlaps_any_range(snps, gtf["exon"]), index=snps.index)

    snps["feature_type"] = "intergenic"
    snps.loc[in_gene, "feature_type"] = "intron"
    snps.loc[in_exon, "feature_type"] = "exon"
    return snps


def apply_region_blacklist_masks(snps, snp_mask, region_bed, blacklist_bed):
    """AND snp_mask with region inclusion and (optional) blacklist exclusion.

    Returns the updated mask and the parsed regions (reused for boundaries).
    """
    regions = read_BED(region_bed)
    region_mask = overlaps_any_range(snps, regions)
    logging.info(f"region filter: {np.sum(region_mask)}/{len(snps)} SNPs passed")
    snp_mask &= region_mask

    if blacklist_bed is not None:
        bl_regions = read_BED(blacklist_bed)
        bl_mask = overlaps_any_range(snps, bl_regions)
        logging.info(
            f"blacklist filter: {np.sum(bl_mask)}/{len(snps)} SNPs in blacklist"
        )
        snp_mask &= ~bl_mask
    return snp_mask, regions


def apply_exon_only_mask(snps, snp_mask, exon_only):
    """Log exonic SNP count and, if exon_only, AND snp_mask with the exon mask."""
    n_exon = int((snps["feature_type"] == "exon").sum())
    logging.info(
        f"#exonic SNPs: {n_exon}/{len(snps)} ({n_exon / max(len(snps), 1):.3%})"
    )
    if exon_only:
        exon_mask = (snps["feature_type"] == "exon").to_numpy()
        logging.info(f"exon filter: {np.sum(exon_mask)}/{len(snps)} SNPs passed")
        snp_mask &= exon_mask
    return snp_mask


def snp_to_region(
    snp_df: pd.DataFrame, region_df: pd.DataFrame, assay_type: str, region_id="BIN_ID"
):
    """Assign SNPs to pre-computed regions and count them per region.

    ``region_df`` must hold 0-based, half-open, non-overlapping intervals. It gains a
    ``#SNPS`` column in place; SNPs outside every region are dropped.
    """
    logging.info(f"#{assay_type}-SNP (raw)={len(snp_df)}")
    snp_df, _ = assign_and_keep(snp_df, region_df, region_id, label=f"{assay_type}: ")
    logging.info(f"#{assay_type}-SNP (remain)={len(snp_df)}")

    counts = snp_df[region_id].value_counts()
    region_df["#SNPS"] = region_df[region_id].map(counts).fillna(0).astype(int)
    return snp_df
