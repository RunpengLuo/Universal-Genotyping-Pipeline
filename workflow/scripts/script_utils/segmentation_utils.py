"""Fixed bins, clustered sums, and the adaptive merge of bins into bbs.

Three groups, in pipeline order:

1. fixed bins    - ``build_fixedwidth_bins`` tiles each segment into ``window_size``
   pieces; one row of the result is a fixed bin, and no bin spans two segments.
2. clustered sums - ``cluster_sum`` and its two wrappers. Every matrix is
   ``(n_features, n_observations)``: binning sums FEATURES (SNPs/bins -> bbs),
   pseudobulking sums OBSERVATIONS (cells -> datasets); both are one sparse one-hot
   multiply.
3. adaptive binning - ``build_adaptive_bins`` walks consecutive fixed bins inside one
   ``seg_id`` and closes a bb once every tumor observation meets its read target. It
   sums the SNP depth per bin through group 2.

The config keys ``min_snp_reads`` / ``min_snp_per_bin`` / ``max_blocksize`` speak of
"bin" in the bb sense; they keep their published names.

GTF-feature annotation and the gene-aware cluster key live in ``feature_utils``.
"""

import logging

import numpy as np
import pandas as pd
import numba

from scipy.sparse import csr_matrix, issparse


##################################################
# 1. fixed bins: tile the segments


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


##################################################
# 2. clustered sums over a count matrix, on either axis


def cluster_sum(X, cluster_ids, n_clusters, axis=0):
    """Sum the features (``axis=0``) or observations (``axis=1``) of *X* within each cluster.

    Args:
        X: ``(n_features, n_observations)`` dense or sparse matrix.
        cluster_ids: Cluster id per feature (``axis=0``) or per observation
            (``axis=1``), in ``[0, n_clusters)``.
        n_clusters: Number of clusters, i.e. the size of the collapsed axis in the output.
        axis: Axis to collapse.

    Returns:
        ``(n_clusters, n_observations)`` for ``axis=0``, ``(n_features, n_clusters)``
        for ``axis=1``; sparse when *X* is sparse.

    Raises:
        ValueError: *cluster_ids* has the wrong length or holds an id outside the range.
    """
    X = X.tocsr() if issparse(X) else np.asarray(X)
    cluster_ids = np.asarray(cluster_ids, dtype=np.int64)
    n = X.shape[axis]
    if cluster_ids.shape[0] != n:
        raise ValueError(
            f"cluster_ids length {cluster_ids.shape[0]} != axis-{axis} size {n}"
        )
    if n and (cluster_ids.min() < 0 or cluster_ids.max() >= n_clusters):
        raise ValueError("cluster_ids out of range")

    onehot = csr_matrix(
        (np.ones(n, dtype=np.int8), (cluster_ids, np.arange(n, dtype=np.int64))),
        shape=(n_clusters, n),
    )
    return onehot @ X if axis == 0 else X @ onehot.T


def sum_features_to_bbs(X, bb_ids, n_bbs):
    """Sum an SNP- or bin-level matrix into ``(n_bbs, n_observations)``."""
    return cluster_sum(X, bb_ids, n_bbs, axis=0)


def sum_observations_to_pseudobulk(mat, cluster_ids, n_clusters):
    """Sum observations into ``(n_features, n_clusters)`` pseudobulks, always dense."""
    out = cluster_sum(mat, cluster_ids, n_clusters, axis=1)
    return out.toarray() if issparse(out) else np.asarray(out)


def dense_observation(mat, i: int):
    """Observation *i* of a dense or sparse matrix, as a 1-D array over features."""
    obs = mat[:, i]
    return obs.toarray().ravel() if issparse(obs) else np.asarray(obs).ravel()


##################################################
# 3. adaptive binning: merge consecutive fixed bins into bbs


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


def build_adaptive_bins(
    bins: pd.DataFrame,
    snps: pd.DataFrame,
    tot_mtx: np.ndarray,
    min_snp_reads,
    min_snp_per_bin: int,
    cluster_cols: list,
    max_blocksize=0,
    gene_aware=False,
):
    """Merge consecutive fixed bins into bbs until the SNP thresholds are met.

    Both *bins* and *snps* gain a ``bb_id`` column IN PLACE, overwritten in full on
    every call, so a caller sweeping several ``min_snp_reads`` reuses the same two
    frames and needs no defensive copy.

    Parameters
    ----------
    bins : pd.DataFrame
        Fixed bins with ``#CHR``, ``START``, ``END``, ``bin_id``, and the clustering
        columns. When ``gene_aware``, must also carry a ``gene_cluster`` column, so a
        bb never splits a gene.
    snps : pd.DataFrame
        SNP DataFrame with ``POS0``, ``#CHR`` and ``bin_id``; the caller assigns the
        fixed bins and drops the SNPs that land in none.
    tot_mtx : (n_snps, M) ndarray
        Per-SNP total read counts over the M TUMOR observations, row-aligned to *snps*;
        the caller slices both.
    min_snp_reads : int or array-like
        Minimum total tumor reads for a bb, per tumor observation. A scalar is
        broadcast to every observation; an array of length M sets one threshold each.
    min_snp_per_bin : int
        Minimum number of SNPs per bb.
    cluster_cols : list of str
        Columns to cluster fixed bins by (e.g. ``["region_id"]``); a bb never spans
        two clusters.

    Returns
    -------
    bbs : pd.DataFrame
        bb definitions with ``#CHR``, ``START``, ``END``, ``#SNPS``, ``BLOCKSIZE``,
        ``bb_id``, and the clustering columns.
    snps : pd.DataFrame
        The *snps* argument itself, now carrying ``bb_id``; rows and order unchanged.
    """
    assert "bin_id" in snps.columns, (
        "snps, no bin_id column; assign the fixed bins first"
    )
    assert len(snps) == tot_mtx.shape[0], (
        f"tot_mtx has {tot_mtx.shape[0]} rows for {len(snps)} SNPs"
    )

    M_tumor = tot_mtx.shape[1]
    min_snp_reads_vec = np.ascontiguousarray(
        np.broadcast_to(np.asarray(min_snp_reads, dtype=np.float64), (M_tumor,))
    )
    logging.info(
        f"build_adaptive_bins: min_snp_reads(per-obs)={min_snp_reads_vec.tolist()}, "
        f"min_snp_per_bin={min_snp_per_bin}, "
        f"max_blocksize={max_blocksize}"
    )

    # 1. Compute per-fixed-bin stats
    B = len(bins)
    snp_bin_ids = snps["bin_id"].to_numpy()

    tot_tumor = (tot_mtx.toarray() if issparse(tot_mtx) else tot_mtx).astype(np.float64)

    bin_nsnps = np.bincount(snp_bin_ids, minlength=B).astype(np.int64)
    bin_reads = np.asarray(
        cluster_sum(tot_tumor, snp_bin_ids, B, axis=0),
        dtype=np.float64,
    )

    # 2. Cluster the fixed bins and run the numba kernel
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

    # 3. Propagate bb_id to the SNPs
    bin_bb_map = bins["bb_id"].to_numpy()
    snps["bb_id"] = bin_bb_map[snps["bin_id"].to_numpy()]

    # 4. Build the bbs frame from the fixed bins
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
