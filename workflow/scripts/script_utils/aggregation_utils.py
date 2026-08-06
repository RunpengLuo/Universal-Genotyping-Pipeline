"""Adaptive binning: merge fixed bins into bbs, and the SNP annotation around it.

A fixed bin is one row of the window BED (a ``window_size`` tile); a bb is the merged
bin that binning emits. ``build_adaptive_bins`` walks consecutive fixed bins inside one
``seg_id`` and closes a bb once every tumor observation meets its read target.

The config keys ``min_snp_reads`` / ``min_snp_per_bin`` / ``max_blocksize`` speak of
"bin" in the bb sense; they keep their published names.
"""

import logging

import numpy as np
import pandas as pd
import numba

from scipy.sparse import issparse

from range_utils import (
    assign_all_features,
    assign_pos_to_range,
    overlaps_any_range,
)
from io_utils import read_BED, read_gtf
from matrix_utils import cluster_sum


def gene_cluster_labels(n_features, ranges):
    """Gene-cluster id per ordered feature (fixed bins or SNPs) so bbs never split a gene.

    Each gene occupies an inclusive index range ``(lo, hi)`` over the feature ordering
    (e.g. the first..last fixed bin holding that gene's SNPs). Every boundary internal to
    a range is made non-cuttable; overlapping ranges (genes that share a feature, i.e.
    adjacent/overlapping genes) merge transitively into one cluster; features in no range
    are singleton clusters (native fixed-bin/SNP granularity). Features sharing a cluster
    id must stay in one bb, so a bb holds whole genes only - never a partial gene.

    Parameters
    ----------
    n_features : int
        Number of ordered features.
    ranges : iterable of (lo, hi)
        Inclusive index ranges, one per gene.

    Returns
    -------
    np.ndarray (int64), length n_features
        Run-length-contiguous cluster id; the boundary between features ``i-1`` and ``i``
        is a cut point iff ``labels[i] != labels[i-1]``.
    """
    if n_features == 0:
        return np.zeros(0, dtype=np.int64)
    blocked = np.zeros(
        n_features - 1, dtype=bool
    )  # blocked[i] = boundary (i, i+1) non-cuttable
    for lo, hi in ranges:
        if hi > lo:
            blocked[lo:hi] = True
    labels = np.empty(n_features, dtype=np.int64)
    labels[0] = 0
    if n_features > 1:
        labels[1:] = np.cumsum(~blocked)
    return labels


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


def assign_and_drop_outside(snps, ref, ref_id, pos_col="POS0", label=""):
    """Assign SNPs to reference ranges, log the misses, and drop them.

    The shared body of ``assign_snps_to_bins`` and ``assign_snps_to_bbs``: assign, report
    the fraction outside every range, drop those SNPs, and cast the id to the reference's
    own dtype.

    Args:
        snps: SNP frame with ``#CHR`` and *pos_col*.
        ref: Reference ranges carrying *ref_id*.
        ref_id: Identifier column to assign.
        pos_col: 0-based position column of *snps*.
        label: Prefix for the log line.

    Returns:
        ``(kept, outside_mask)``: the assigned SNPs reindexed from 0, and the
        boolean mask of dropped SNPs over the INPUT frame.
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


def assign_snps_to_bins(snps, bins, tot_mtx):
    """Assign SNPs to fixed bins once, dropping those outside every bin.

    Hoisted out of ``build_adaptive_bins`` so a caller sweeping binning parameters pays
    for the assignment once instead of per sweep point. The returned frame carries
    ``_orig_df_idx`` (the SNP's position in the input frame, indexing ``tot_mtx``) and an
    int64 ``bin_id``; ``build_adaptive_bins`` skips step 1 when it sees one.

    Args:
        snps: SNP frame with ``#CHR`` and ``POS0``.
        bins: Fixed-bin frame with ``#CHR``, ``START``, ``END``, ``bin_id``.
        tot_mtx: Per-SNP total counts, for the off-bin depth log only.

    Returns:
        The SNPs inside a fixed bin, reindexed from 0.
    """
    snps["_orig_df_idx"] = np.arange(len(snps))
    orig_df_idx = snps["_orig_df_idx"].to_numpy()
    kept, outside = assign_and_drop_outside(snps, bins, "bin_id")
    if outside.any():
        off_depth = tot_mtx[orig_df_idx[outside]].sum(axis=1)
        logging.info(
            f"off-target SNP depth: "
            f"min={off_depth.min()}, max={off_depth.max()}, "
            f"mean={off_depth.mean():.1f}, median={np.median(off_depth):.1f}"
        )
    kept["bin_id"] = kept["bin_id"].astype(np.int64)
    return kept


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
        ``gene_cluster_labels``) so bbs never split a gene.
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

    # 1. Assign SNPs to fixed bins, unless the caller already did (assign_snps_to_bins)
    if "bin_id" not in snps.columns:
        snps = assign_snps_to_bins(snps, bins, tot_mtx)

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

    # PS column if present
    if "PS" in snps.columns:
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

    Returns the updated mask and the parsed regions (reused for the SNP ranges).
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


def assign_snps_to_bbs(
    snp_df: pd.DataFrame, bb_df: pd.DataFrame, assay_type: str, id_col="bb_id"
):
    """Assign SNPs to pre-computed bbs and count them per bb.

    ``bb_df`` must hold 0-based, half-open, non-overlapping ranges. It gains a ``#SNPS``
    column in place; SNPs outside every bb are dropped.
    """
    logging.info(f"#{assay_type}-SNP (raw)={len(snp_df)}")
    snp_df, _ = assign_and_drop_outside(snp_df, bb_df, id_col, label=f"{assay_type}: ")
    logging.info(f"#{assay_type}-SNP (remain)={len(snp_df)}")

    counts = snp_df[id_col].value_counts()
    bb_df["#SNPS"] = bb_df[id_col].map(counts).fillna(0).astype(int)
    return snp_df
