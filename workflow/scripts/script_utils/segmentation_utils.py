"""Fixed bins, clustered sums, and the adaptive merge of bins into bbs.

Last update: 2026-08-11

Functions:
- build_fixedwidth_bins: tile segments, so no bin spans two
- cluster_sum: collapse one matrix axis by a cluster-id array
- sum_features_to_bbs: the feature-axis wrapper, SNPs or bins into bbs
- sum_observations_to_pseudobulk: the observation-axis wrapper, cells into datasets
- dense_observation: one matrix column as a 1-D array
- build_adaptive_bins: merge bins until min_snp_reads and min_total_reads are met
References:
- docs/DEVELOPER.md: bin is the fixed tile, bb the merged one
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
    bin_total_reads,
    min_total_reads,
    require_snp,
    cluster_ids,
):
    """Greedy adaptive merge of consecutive fixed bins into bbs.

    A bb closes at the first bin where both criteria hold::

        snp_ok   = acc[j] >= min_snp_reads_vec[j] for every tumor j, and
                   acc_n >= min_snp_per_bin
        total_ok = acc_tot[j] >= min_total_reads for every observation j

        close when (snp_ok or not require_snp) and total_ok

    There is no span cap: a cap can only fire by cutting a bin that has not met the read
    criteria, which contradicts them. The only bb that may fall short is a cluster's
    trailing run, which merges back into the previous bb; when the whole cluster is one
    run it stays one bb, as HATCHet2 does per chromosome arm (``adaptive_bins_arm``,
    doi:10.1038/s41587-020-0661-6).

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
    bin_total_reads : (B, M_all) contiguous float64
        Per-fixed-bin total reads over EVERY observation, normal included. Read starts,
        not depth: an integer count with no missing entries.
    min_total_reads : float
        Minimum total reads per bb per observation; 0 disables the criterion.
    require_snp : bool
        False drops the SNP criterion, leaving *min_total_reads* alone to close the bb.
        Set inside a clonal-LOH region, where a germline het is masked and no SNP count
        means anything.
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
    M_all = bin_total_reads.shape[1]
    bb_ids = np.zeros(B, dtype=np.int64)
    if B == 0:
        return bb_ids, 0

    bb_id = 0
    prev_start = 0
    acc = bin_reads[0].copy()
    acc_tot = bin_total_reads[0].copy()
    acc_n = bin_nsnps[0]

    for i in range(1, B):
        snp_ok = acc_n >= min_snp_per_bin
        if snp_ok:
            for j in range(M):
                if acc[j] < min_snp_reads_vec[j]:
                    snp_ok = False
                    break
        total_ok = True
        for j in range(M_all):
            if acc_tot[j] < min_total_reads:
                total_ok = False
                break
        cluster_boundary = cluster_ids[i] != cluster_ids[i - 1]
        if (snp_ok or not require_snp) and total_ok and cluster_boundary:
            bb_ids[prev_start:i] = bb_id
            bb_id += 1
            prev_start = i
            acc = bin_reads[i].copy()
            acc_tot = bin_total_reads[i].copy()
            acc_n = bin_nsnps[i]
        else:
            for j in range(M):
                acc[j] += bin_reads[i, j]
            for j in range(M_all):
                acc_tot[j] += bin_total_reads[i, j]
            acc_n += bin_nsnps[i]

    # last run: its own bb when it meets the same test, else merged into the previous one
    last_snp_ok = acc_n >= min_snp_per_bin
    if last_snp_ok:
        for j in range(M):
            if acc[j] < min_snp_reads_vec[j]:
                last_snp_ok = False
                break
    last_total_ok = True
    for j in range(M_all):
        if acc_tot[j] < min_total_reads:
            last_total_ok = False
            break
    if (last_snp_ok or not require_snp) and last_total_ok and prev_start > 0:
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
    gene_aware=False,
    bin_total_reads=None,
    min_total_reads=0,
    is_loh_col=None,
):
    """Merge consecutive fixed bins into bbs until the read thresholds are met.

    A bb closes when its SNP reads meet ``min_snp_reads`` in every tumor AND its total
    reads meet ``min_total_reads`` in every observation; see :func:`_merge_bins_to_bbs`
    for the exact rule. Inside a clonal-LOH region the SNP criterion is dropped, since
    the germline hets there are masked.

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
    gene_aware : bool
        Close a bb only at a ``gene_cluster`` boundary.
    bin_total_reads : (B, M_all) ndarray or None
        Per-fixed-bin read starts over every observation, row-aligned to *bins*. None
        leaves the total-read criterion off, whatever *min_total_reads* says.
    min_total_reads : float
        Minimum read starts per bb per observation; 0 disables the criterion.
    is_loh_col : str or None
        Boolean column of *bins* marking clonal LOH. Must be constant within a cluster,
        which holds when the caller puts its LOH id in *cluster_cols*.

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

    # 1. Compute per-fixed-bin stats
    B = len(bins)
    snp_bin_ids = snps["bin_id"].to_numpy()

    tot_tumor = (tot_mtx.toarray() if issparse(tot_mtx) else tot_mtx).astype(np.float64)

    bin_nsnps = np.bincount(snp_bin_ids, minlength=B).astype(np.int64)
    bin_reads = np.asarray(
        cluster_sum(tot_tumor, snp_bin_ids, B, axis=0),
        dtype=np.float64,
    )
    if bin_total_reads is None:
        min_total_reads = 0.0
        bin_totals = np.zeros((B, 1), dtype=np.float64)
    else:
        bin_totals = np.ascontiguousarray(
            (
                bin_total_reads.toarray()
                if issparse(bin_total_reads)
                else np.asarray(bin_total_reads)
            ),
            dtype=np.float64,
        )
        assert bin_totals.shape[0] == B, (
            f"bin_total_reads has {bin_totals.shape[0]} rows for {B} fixed bins"
        )
    min_total_reads = float(min_total_reads)
    logging.info(
        f"build_adaptive_bins: min_snp_reads(per-obs)={min_snp_reads_vec.tolist()}, "
        f"min_snp_per_bin={min_snp_per_bin}, min_total_reads={min_total_reads:g}"
    )

    # 2. Cluster the fixed bins and run the numba kernel
    bb_id = 0
    bins["bb_id"] = 0
    # groupby drops null keys, which would leave those bins on the initialized bb_id 0
    assert bins[cluster_cols].notna().all().all(), (
        f"fixed bins, null cluster key in {cluster_cols}"
    )
    is_loh = (
        bins[is_loh_col].to_numpy(dtype=bool)
        if is_loh_col is not None
        else np.zeros(B, dtype=bool)
    )
    bin_clusters = bins.groupby(by=cluster_cols, sort=False)
    logging.info(f"#fixed-bin clusters={len(bin_clusters)}, keys: {cluster_cols}")

    n_loh_clusters = 0
    for _, cluster_bins in bin_clusters:
        idxs = cluster_bins.index.to_numpy()
        c_reads = np.ascontiguousarray(bin_reads[idxs])
        c_nsnps = np.ascontiguousarray(bin_nsnps[idxs])
        c_totals = np.ascontiguousarray(bin_totals[idxs])
        c_loh = is_loh[idxs]
        assert c_loh.all() == c_loh.any(), (
            f"cluster mixes LOH and non-LOH bins; put the LOH id in {cluster_cols}"
        )
        require_snp = not bool(c_loh[0])
        n_loh_clusters += 0 if require_snp else 1
        if gene_aware:
            c_gene = cluster_bins["gene_cluster"].to_numpy()
        else:
            c_gene = np.arange(len(idxs), dtype=np.int64)

        local_bb_ids, n_bbs = _merge_bins_to_bbs(
            c_reads,
            c_nsnps,
            min_snp_reads_vec,
            min_snp_per_bin,
            c_totals,
            min_total_reads,
            require_snp,
            c_gene,
        )

        local_bb_ids += bb_id
        bins.loc[idxs, "bb_id"] = local_bb_ids
        bb_id += max(n_bbs, 1)

    if is_loh_col is not None:
        logging.info(
            f"clusters binned without the SNP criterion (clonal LOH): "
            f"{n_loh_clusters}/{len(bin_clusters)}"
        )

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

    # bbs that ran out of cluster before meeting a threshold: HATCHet2's one-bin-per-arm
    # case. A spike here means the thresholds are too high for the coverage.
    bb_of_bin = bins["bb_id"].to_numpy()
    loh_bb = np.zeros(num_bbs, dtype=bool)
    loh_bb[bb_of_bin[is_loh]] = True
    bb_reads = np.asarray(cluster_sum(bin_reads, bb_of_bin, num_bbs, axis=0))
    short_snp = (bb_reads < min_snp_reads_vec).any(axis=1) & ~loh_bb
    logging.info(
        f"bbs below min_snp_reads: {int(short_snp.sum())}/{num_bbs} "
        f"(cluster ended first; LOH bbs excluded)"
    )
    if min_total_reads > 0:
        bb_totals = np.asarray(cluster_sum(bin_totals, bb_of_bin, num_bbs, axis=0))
        short_tot = (bb_totals < min_total_reads).any(axis=1)
        logging.info(
            f"bbs below min_total_reads: {int(short_tot.sum())}/{num_bbs} "
            f"(cluster ended first)"
        )

    return bbs, snps
