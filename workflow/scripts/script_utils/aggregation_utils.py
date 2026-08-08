"""Adaptive binning: merge fixed bins into bbs, and the SNP filtering around it.

A fixed bin is one row of the window BED (a ``window_size`` tile); a bb is the merged
bin that binning emits. ``build_adaptive_bins`` walks consecutive fixed bins inside one
``seg_id`` and closes a bb once every tumor observation meets its read target.

The config keys ``min_snp_reads`` / ``min_snp_per_bin`` / ``max_blocksize`` speak of
"bin" in the bb sense; they keep their published names.

GTF-feature annotation and the gene-aware cluster key live in ``feature_utils``.
"""

import logging

import numpy as np
import pandas as pd
import numba

from scipy.sparse import issparse

from range_utils import assign_pos_to_range
from matrix_utils import cluster_sum


def build_fixedwidth_bins(segments, bin_size, chroms=None):
    """Tile every segment into fixed-width bins, so no bin spans two segments.

    Each segment ``[START, END)`` is cut into *bin_size* pieces; the trailing remainder
    becomes one shorter bin, unless it is under half a bin, in which case it is absorbed
    by its predecessor. Every column of *segments* beyond the coordinates (``region_id``,
    ``seg_id``, ...) is carried onto the bins it produced, so the ids need no separate
    assignment pass.

    Fully vectorized: bin counts come from one ``divmod`` over the segments and the rows
    from one ``repeat``, so cost is O(#bins) in numpy rather than a Python loop per bin.

    Args:
        segments: Segments with ``#CHR``, ``START``, ``END``, plus any id columns.
        bin_size: Bin width in bp.
        chroms: Keep only segments on these contigs; ``None`` keeps every segment.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END`` and the carried columns, in segment
        order.
    """
    if chroms is not None:
        segments = segments[segments["#CHR"].isin(chroms)]
    starts = segments["START"].to_numpy(dtype=np.int64)
    ends = segments["END"].to_numpy(dtype=np.int64)

    n_full, rem = np.divmod(ends - starts, bin_size)
    # a trailing piece under half a bin joins its predecessor instead of standing alone
    stands_alone = (rem > 0) & ~((n_full >= 1) & (rem < bin_size // 2))
    n_bins = n_full + stands_alone

    seg_idx = np.repeat(np.arange(len(segments)), n_bins)
    # index of each bin within its own segment
    k = np.arange(len(seg_idx)) - np.repeat(np.cumsum(n_bins) - n_bins, n_bins)

    bin_starts = starts[seg_idx] + k * bin_size
    # the last bin of a segment always closes on the segment's own END, whether it is
    # short (remainder) or long (absorbed remainder)
    bin_ends = np.where(k == n_bins[seg_idx] - 1, ends[seg_idx], bin_starts + bin_size)

    out = {"#CHR": segments["#CHR"].to_numpy()[seg_idx]}
    out["START"], out["END"] = bin_starts, bin_ends
    for col in segments.columns:
        if col not in ("#CHR", "START", "END"):
            out[col] = segments[col].to_numpy()[seg_idx]
    return pd.DataFrame(out)


def stamp_bin_label(bin_df, snps_binned, col, default):
    """Carry a per-SNP label onto the fixed bins as a never-null cluster key.

    Each bin takes the modal value of *col* over its SNPs, then bins with no SNP inherit
    from their neighbours (forward, then backward), and a frame with no SNP at all falls
    back to *default*. The fill must cover both ends: ``build_adaptive_bins`` groups the
    bins by ``cluster_cols`` and pandas drops null keys, so a null-keyed bin would never
    enter the merge loop and would keep the initialized ``bb_id`` of 0.

    Args:
        bin_df: Fixed bins with ``bin_id``. Modified in place.
        snps_binned: SNPs carrying ``bin_id`` and *col*.
        col: Label column to carry over.
        default: Value for bins when no SNP in the frame has one.

    Returns:
        *bin_df* with *col* added.
    """
    per_bin = snps_binned.groupby("bin_id")[col].agg(lambda x: x.mode().iloc[0])
    bin_df[col] = bin_df["bin_id"].map(per_bin)
    if bin_df[col].isna().any():
        bin_df[col] = bin_df[col].ffill().bfill().fillna(default)
    return bin_df


@numba.njit
def _merge_bins_to_bbs(
    bin_reads,
    bin_nsnps,
    min_snp_reads_vec,
    min_snp_per_bin,
    bin_starts,
    bin_ends,
    max_blocksize,
    cluster_ids,
):
    """Greedy adaptive merge of consecutive fixed bins into bbs.

    A bb can close only when its accumulated reads meet the per-observation threshold
    (``acc[j] >= min_snp_reads_vec[j]`` for every observation j). ``max_blocksize`` is a
    SECONDARY guard: it may force a cut only once that read threshold is already
    satisfied, so no bb is ever emitted below threshold and every bb is non-empty in
    every observation (e.g. every WES observation carries reads).

    Parameters
    ----------
    bin_reads : (B, M) contiguous float64
        Per-fixed-bin total tumor reads (summed from the SNPs in each bin).
    bin_nsnps : (B,) int64
        Number of SNPs per fixed bin.
    min_snp_reads_vec : (M,) float64
        Per-observation minimum SNP reads per bb.
    min_snp_per_bin : int
        Minimum number of SNPs per bb.
    bin_starts : (B,) int64
        START coordinate per fixed bin.
    bin_ends : (B,) int64
        END coordinate per fixed bin.
    max_blocksize : int
        Span cap (bp); once the read threshold is met, a bb over this span is cut even
        if it holds fewer than ``min_snp_per_bin`` SNPs. Never cuts below the read
        threshold. Set to 0 to disable.
    cluster_ids : (B,) int64
        Gene-cluster id per fixed bin; a bb may only close at a cluster boundary
        (``cluster_ids[i] != cluster_ids[i-1]``), so a gene is never split across bbs.
        Pass ``np.arange(B)`` for no constraint.

    Returns
    -------
    bb_ids : (B,) int64
        Relative bb ID for each fixed bin within this cluster.
    n_bbs : int
        Number of bbs created (before last-run adjustment).
    """
    B, M = bin_reads.shape
    bb_ids = np.zeros(B, dtype=np.int64)
    if B == 0:
        return bb_ids, 0

    bb_id = 0
    prev_start = 0
    acc = bin_reads[0].copy()
    acc_n = bin_nsnps[0]

    for i in range(1, B):
        meets_reads = True
        for j in range(M):
            if acc[j] < min_snp_reads_vec[j]:
                meets_reads = False
                break
        span = bin_ends[i - 1] - bin_starts[prev_start]
        exceeds_size = max_blocksize > 0 and span >= max_blocksize
        cluster_boundary = cluster_ids[i] != cluster_ids[i - 1]
        if (
            meets_reads
            and (acc_n >= min_snp_per_bin or exceeds_size)
            and cluster_boundary
        ):
            bb_ids[prev_start:i] = bb_id
            bb_id += 1
            prev_start = i
            acc = bin_reads[i].copy()
            acc_n = bin_nsnps[i]
        else:
            for j in range(M):
                acc[j] += bin_reads[i, j]
            acc_n += bin_nsnps[i]

    # last run: keep as its own bb if it meets the read threshold and isn't the only run
    last_meets_reads = True
    for j in range(M):
        if acc[j] < min_snp_reads_vec[j]:
            last_meets_reads = False
            break
    last_span = bin_ends[B - 1] - bin_starts[prev_start]
    last_exceeds_size = max_blocksize > 0 and last_span >= max_blocksize
    if (
        last_meets_reads
        and (acc_n >= min_snp_per_bin or last_exceeds_size)
        and prev_start > 0
    ):
        bb_ids[prev_start:] = bb_id
        bb_id += 1
    else:
        merge_id = bb_id - 1 if bb_id > 0 else 0
        bb_ids[prev_start:] = merge_id

    return bb_ids, bb_id


def log_off_range_depth(na_idx, tot_mtx, label=""):
    """Log the depth carried by the positions an assignment dropped.

    *na_idx* is the second return of any ``assign_*`` and indexes *tot_mtx* rows
    directly, so this reports how much signal falls outside the reference grid.

    Args:
        na_idx: Positional indices of the unassigned rows.
        tot_mtx: Per-position total counts, rows aligned to the pre-assignment frame.
        label: Prefix for the log line.
    """
    if not len(na_idx):
        return
    off_depth = tot_mtx[na_idx].sum(axis=1)
    logging.info(
        f"{label}off-range depth over {len(na_idx)} dropped positions: "
        f"min={off_depth.min()}, max={off_depth.max()}, "
        f"mean={off_depth.mean():.1f}, median={np.median(off_depth):.1f}"
    )


def build_adaptive_bins(
    bins: pd.DataFrame,
    snps: pd.DataFrame,
    tot_mtx: np.ndarray,
    min_snp_reads,
    min_snp_per_bin: int,
    cluster_cols: list,
    tumor_sidx=0,
    max_blocksize=0,
    gene_aware=False,
):
    """Merge consecutive fixed bins into bbs until the SNP thresholds are met.

    Parameters
    ----------
    bins : pd.DataFrame
        Fixed bins with ``#CHR``, ``START``, ``END``, ``bin_id``, and the clustering
        columns. When ``gene_aware``, must also carry a ``gene_cluster`` column (see
        ``feature_utils.stamp_gene_clusters``) so bbs never split a gene.
    snps : pd.DataFrame
        SNP DataFrame with ``POS0`` and ``#CHR`` columns.
    tot_mtx : (n_snps, M) ndarray
        Per-SNP total read counts over M observations.
    min_snp_reads : int or array-like
        Minimum total tumor reads for a bb, per tumor observation. A scalar is
        broadcast to every tumor observation; an array of length ``M - tumor_sidx``
        sets a per-observation threshold.
    min_snp_per_bin : int
        Minimum number of SNPs per bb.
    cluster_cols : list of str
        Columns to cluster fixed bins by (e.g. ``["region_id"]``); a bb never spans
        two clusters.
    tumor_sidx : int
        Index of the first tumor observation.

    Returns
    -------
    bbs : pd.DataFrame
        bb definitions with ``#CHR``, ``START``, ``END``, ``#SNPS``, ``BLOCKSIZE``,
        ``bb_id``, and the clustering columns.
    snps : pd.DataFrame
        Input SNPs with ``bb_id`` and ``bin_id`` columns added. SNPs not falling in
        any fixed bin are dropped.
    """

    M_tumor = tot_mtx.shape[1] - tumor_sidx
    min_snp_reads_vec = np.ascontiguousarray(
        np.broadcast_to(np.asarray(min_snp_reads, dtype=np.float64), (M_tumor,))
    )
    logging.info(
        f"build_adaptive_bins: min_snp_reads(per-obs)={min_snp_reads_vec.tolist()}, "
        f"min_snp_per_bin={min_snp_per_bin}, "
        f"max_blocksize={max_blocksize}"
    )

    # 1. Assign SNPs to fixed bins, unless the caller already did
    if "bin_id" not in snps.columns:
        snps["_orig_df_idx"] = np.arange(len(snps))
        snps, na_idx = assign_pos_to_range(snps, bins, ref_id="bin_id", dropna=True)
        log_off_range_depth(na_idx, tot_mtx)

    # 2. Compute per-fixed-bin stats
    B = len(bins)
    snp_bin_ids = snps["bin_id"].to_numpy()
    snp_orig_df_idx = snps["_orig_df_idx"].to_numpy()

    if issparse(tot_mtx):
        tot_tumor = tot_mtx[:, tumor_sidx:].toarray().astype(np.float64)
    else:
        tot_tumor = tot_mtx[:, tumor_sidx:].astype(np.float64)

    bin_nsnps = np.bincount(snp_bin_ids, minlength=B).astype(np.int64)
    bin_reads = np.asarray(
        cluster_sum(tot_tumor[snp_orig_df_idx], snp_bin_ids, B, axis=0),
        dtype=np.float64,
    )

    # 3. Cluster the fixed bins and run the numba kernel
    bb_id = 0
    bins["bb_id"] = 0
    all_bin_starts = bins["START"].to_numpy(dtype=np.int64)
    all_bin_ends = bins["END"].to_numpy(dtype=np.int64)
    # groupby drops null keys, which would leave those bins on the initialized bb_id 0
    assert bins[cluster_cols].notna().all().all(), (
        f"fixed bins, null cluster key in {cluster_cols}"
    )
    bin_clusters = bins.groupby(by=cluster_cols, sort=False)
    logging.info(f"#fixed-bin clusters={len(bin_clusters)}, keys: {cluster_cols}")

    for _, cluster_bins in bin_clusters:
        idxs = cluster_bins.index.to_numpy()
        c_reads = np.ascontiguousarray(bin_reads[idxs])
        c_nsnps = np.ascontiguousarray(bin_nsnps[idxs])
        c_starts = np.ascontiguousarray(all_bin_starts[idxs])
        c_ends = np.ascontiguousarray(all_bin_ends[idxs])
        if gene_aware:
            c_gene = cluster_bins["gene_cluster"].to_numpy()
        else:
            c_gene = np.arange(len(idxs), dtype=np.int64)

        local_bb_ids, n_bbs = _merge_bins_to_bbs(
            c_reads,
            c_nsnps,
            min_snp_reads_vec,
            min_snp_per_bin,
            c_starts,
            c_ends,
            max_blocksize,
            c_gene,
        )

        local_bb_ids += bb_id
        bins.loc[idxs, "bb_id"] = local_bb_ids
        bb_id += max(n_bbs, 1)

    # 4. Propagate bb_id to the SNPs
    bin_bb_map = bins["bb_id"].to_numpy()
    snps["bb_id"] = bin_bb_map[snps["bin_id"].to_numpy()]

    # 5. Build the bbs frame from the fixed bins
    pos_dict = {
        "#CHR": ("#CHR", "first"),
        "START": ("START", "min"),
        "END": ("END", "max"),
    }
    for col in cluster_cols:
        pos_dict[col] = (col, "first")

    # drop the index name: it would collide with the bb_id column added below, and
    # pandas rejects `join(on="bb_id")` when the name is both an index level and a column
    bbs = bins.groupby("bb_id", sort=True).agg(**pos_dict).rename_axis(None)

    # SNP counts per bb
    snp_counts = snps.groupby("bb_id").size()
    bbs["#SNPS"] = bbs.index.map(snp_counts).fillna(0).astype(int)
    bbs["BLOCKSIZE"] = bbs["END"] - bbs["START"]
    bbs["bb_id"] = bbs.index

    # PS column if present; a cluster_cols PS is already carried by pos_dict, and taking
    # it from the SNPs instead would be NaN for a bb that holds none
    if "PS" in snps.columns and "PS" not in bbs.columns:
        ps = snps.groupby("bb_id")["PS"].first()
        bbs["PS"] = bbs.index.map(ps)

    num_bbs = len(bbs)
    bb_nsnps = bbs["#SNPS"].to_numpy()
    bb_spans = bbs["BLOCKSIZE"].to_numpy()
    logging.info("build_adaptive_bins summary")
    logging.info(f"#SNPs={len(snps)}, #fixed bins={B}, #bbs={num_bbs}")
    logging.info(
        "snps per bb: min=%.0f  median=%.0f  max=%.0f",
        float(bb_nsnps.min()) if num_bbs > 0 else 0,
        float(np.median(bb_nsnps)) if num_bbs > 0 else 0,
        float(bb_nsnps.max()) if num_bbs > 0 else 0,
    )
    logging.info(
        "span per bb: min=%.0f  median=%.0f  max=%.0f",
        float(bb_spans.min()) if num_bbs > 0 else 0,
        float(np.median(bb_spans)) if num_bbs > 0 else 0,
        float(bb_spans.max()) if num_bbs > 0 else 0,
    )

    return bbs, snps
