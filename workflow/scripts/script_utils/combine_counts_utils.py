"""Observation bookkeeping and bulk depth/RDR summaries for combine_counts.

Last update: 2026-08-13

Functions:
- observation_cluster_ids: map each matrix column to its roster row
- tumor_observation_indices: the columns min_snp_reads is thresholded against
- summarize_observation_counts: per-dataset unit-level QC rows for a single-cell assay
- summarize_read_depth_bb: length-weighted fixed-bin depth aggregated per bb
- summarize_rdr_bb: per-bb RDR, matched-normal or median-centered
"""

import logging

import numpy as np
import pandas as pd

from utils import log_ratios


##################################################
# observations: which matrix column is which replicate/assay


OBS_KEY = ["dataset_id", "assay_type"]


def observation_cluster_ids(cells: pd.DataFrame, roster: pd.DataFrame):
    """Cluster id per matrix column: its row index in the *roster*.

    A multiome pair shares one ``dataset_id``, so the cluster key is the
    ``(dataset_id, assay_type)`` pair, not the dataset alone.

    Args:
        cells: One row per matrix column, from ``read_barcodes_by_dataset``.
        roster: ``sample_ids.tsv``, one row per pseudobulk observation; its row order
            defines the codes.

    Returns:
        int64 index into *roster*, one per column of *cells*.
    """
    keys = pd.MultiIndex.from_frame(roster[OBS_KEY].astype(str))
    assert keys.is_unique, "sample_ids.tsv, duplicate (dataset_id, assay_type)"
    codes = keys.get_indexer(pd.MultiIndex.from_frame(cells[OBS_KEY].astype(str)))
    assert (codes >= 0).all(), (
        "barcodes, (dataset_id, assay_type) absent from sample_ids.tsv"
    )
    return codes.astype(np.int64)


def tumor_observation_indices(roster: pd.DataFrame):
    """Matrix columns of the tumor observations, in column order.

    ``min_snp_reads`` is a per-tumor read target, so binning thresholds these columns
    only: a normal serves as an RDR base or a reference and its BAF is uninformative.
    Both modes select the binning columns through here, so bulk and single-cell agree.

    Args:
        roster: ``sample_ids.tsv``, one row per observation, carrying ``sample_type``.

    Returns:
        Indices into the observation axis; empty when no observation is a tumor, which
        leaves ``build_adaptive_bins`` with no read threshold at all.
    """
    idx = np.flatnonzero(roster["sample_type"].to_numpy() == "tumor").tolist()
    if not idx:
        logging.warning(
            f"no tumor observation among {len(roster)} rows; min_snp_reads has no column "
            "to threshold and binning falls back to min_snp_per_bin alone"
        )
    return idx


##################################################
# single-cell unit-level QC


def summarize_observation_counts(
    sample_df, sample_id, observation_ids, snp_mtx, unit_mtx
):
    """Per-dataset unit-level QC over one assay's cells/spots.

    Both matrices are ``(n_features, n_observations)`` on the same observation axis:
    *snp_mtx* the SNP-covered counts over the shared SNP grid and *unit_mtx* the native
    counts over the assay's own unit - genes for RNA, windows for scATAC.

    ``#SNP`` and ``#feature`` are the rows this dataset actually detects, not the grid
    size: a locus or feature with at least one read in at least one of its cells. The
    grid size is the same for every row and is logged instead. The two sparsities are
    the zero fraction over those detected rows only, so an undetected locus does not
    count twice - once by lowering ``#SNP`` and again by inflating the sparsity.
    ``mean_*_reads`` and ``median_*_reads`` are over that dataset's per-cell/spot
    totals.

    Args:
        sample_df: ``sample_ids.tsv`` for this assay, one row per dataset.
        sample_id: Patient/sample id, stamped on every row.
        observation_ids: Length-n_observations index of each column into *sample_df*.
        snp_mtx: (n_snps, n_observations) SNP-covered counts, sparse or dense.
        unit_mtx: (n_features, n_observations) native counts, sparse or dense.

    Returns:
        ``(stats_df, snp_totals, unit_totals)``. *stats_df* has one row per *sample_df*
        row; the totals are per-observation column sums of each matrix, float64.
    """
    snp_totals = _observation_totals(snp_mtx)
    unit_totals = _observation_totals(unit_mtx)
    observation_ids = np.asarray(observation_ids)
    logging.info(
        f"unit-level QC grid: {snp_mtx.shape[0]} SNPs, {unit_mtx.shape[0]} features "
        f"over {len(observation_ids)} observations"
    )
    rows = []
    for i, row in enumerate(sample_df.itertuples(index=False)):
        cols = np.flatnonzero(observation_ids == i)
        n_snp, snp_sparsity = _detected_and_sparsity(snp_mtx, cols)
        n_feature, feature_sparsity = _detected_and_sparsity(unit_mtx, cols)
        rows.append(
            {
                "sample_id": sample_id,
                "dataset_id": row.dataset_id,
                "sample_type": row.sample_type,
                "assay_type": row.assay_type,
                "#observation": len(cols),
                "#SNP": n_snp,
                "#feature": n_feature,
                "snp_sparsity": snp_sparsity,
                "feature_sparsity": feature_sparsity,
                "mean_snp_reads": _stat(np.mean, snp_totals, cols),
                "median_snp_reads": _stat(np.median, snp_totals, cols),
                "mean_feature_reads": _stat(np.mean, unit_totals, cols),
                "median_feature_reads": _stat(np.median, unit_totals, cols),
            }
        )
    return pd.DataFrame(rows), snp_totals, unit_totals


def _observation_totals(mtx):
    """Per-observation column sums as a float64 1-D array."""
    total = mtx.sum(axis=0)
    return np.asarray(total, dtype=np.float64).ravel()


def _detected_and_sparsity(mtx, cols):
    """Rows of ``mtx[:, cols]`` carrying a read, and the zero fraction among them.

    Rows this dataset never covers are excluded from both, so the sparsity measures how
    thinly the detected rows are populated across its observations.
    """
    if len(cols) == 0:
        return 0, 1.0
    sub = mtx[:, cols]
    per_row = sub.getnnz(axis=1) if hasattr(sub, "getnnz") else np.count_nonzero(sub, 1)
    n_detected = int(np.count_nonzero(per_row))
    if n_detected == 0:
        return 0, 1.0
    return n_detected, float(1.0 - int(per_row.sum()) / (n_detected * len(cols)))


def _stat(fn, totals, cols):
    """*fn* over this dataset's per-observation totals; 0.0 when it owns no column."""
    return float(fn(totals[cols])) if len(cols) else 0.0


##################################################
# bulk read depth and RDR


def summarize_read_depth_bb(bin_df, dp_corrected, num_bbs, dataset_ids=None):
    """Length-weighted aggregation of corrected fixed-bin depth into bbs.

    *bin_df* is the shared fixed-bin grid carrying the ``bb_id`` that
    ``build_adaptive_bins`` stamped, row-aligned to *dp_corrected*. NaN is masked per
    column: a window the correction could not resolve for one dataset still counts for
    every other dataset.

    Args:
        bin_df: Fixed bins with ``START``, ``END`` and ``bb_id``.
        dp_corrected: ``(n_bins, n_datasets)`` corrected depth, NaN where undefined.
        num_bbs: Number of bbs.
        dataset_ids: Column labels for the log line; falls back to the column index.

    Returns:
        ``(bb_dp, bb_bases)``: per-bb mean depth and per-bb total aligned bases, both
        ``(num_bbs, n_datasets)``; *bb_dp* is NaN for a bb with no finite window.
    """
    bb_ids = bin_df["bb_id"].to_numpy(np.int64)
    bin_lengths = (bin_df["END"] - bin_df["START"]).to_numpy(dtype=np.float64)
    num_datasets = dp_corrected.shape[1]
    n_bins = dp_corrected.shape[0]
    bb_dp = np.full((num_bbs, num_datasets), np.nan, dtype=np.float32)
    bb_bases = np.zeros((num_bbs, num_datasets), dtype=np.float64)
    for s in range(num_datasets):
        label = dataset_ids[s] if dataset_ids is not None else f"column {s}"
        finite = np.isfinite(dp_corrected[:, s])
        n_nan = len(finite) - int(finite.sum())
        log_ratios(
            f"{label} NaN fixed bins masked",
            n_nan,
            n_bins,
            bin_lengths[~finite],
            "bin",
            prefix="  ",
        )
        ids, lengths = bb_ids[finite], bin_lengths[finite]
        weighted_sums = np.bincount(
            ids, weights=dp_corrected[finite, s] * lengths, minlength=num_bbs
        )
        total_len_per_bb = np.bincount(ids, weights=lengths, minlength=num_bbs)
        bb_bases[:, s] = weighted_sums
        with np.errstate(invalid="ignore", divide="ignore"):
            bb_dp[:, s] = weighted_sums / total_len_per_bb
    nan_bb = np.isnan(bb_dp).any(axis=1)
    bb_lengths = np.bincount(bb_ids, weights=bin_lengths, minlength=num_bbs)
    log_ratios(
        "bb depth NaN (no finite window in some column)",
        int(nan_bb.sum()),
        num_bbs,
        bb_lengths[nan_bb],
        "bb",
        prefix="  ",
    )
    return bb_dp, bb_bases


def summarize_rdr_bb(
    bin_df,
    dp_corrected,
    bb_dp,
    tumor_dataset_indices,
    get_rdr_base_dataset_id,
    dataset_ids,
):
    """Per-bb RDR for every tumor observation.

    Each tumor with an RDR base observation (``{tumor index: base index}``) is normalized
    by that base, library-size corrected; a tumor without a base is median-centered. The base
    may be any observation (e.g. a different assay/platform), so library sizes are
    computed globally per observation over the shared fixed-bin grid. Returns a
    ``(num_bbs, len(tumor_dataset_indices))`` array aligned to *tumor_dataset_indices*.
    A bb whose base has zero depth yields NaN, not an infinite ratio. High RDR is left
    alone: a focal amplification is signal, and unmappable sequence is already masked.
    """
    num_bbs = bb_dp.shape[0]
    bb_rdr = np.full((num_bbs, len(tumor_dataset_indices)), np.nan, dtype=np.float32)
    rdr_pos = {o: i for i, o in enumerate(tumor_dataset_indices)}

    # global per-observation total aligned bases for library-size correction
    bin_sizes = (bin_df["END"] - bin_df["START"]).to_numpy(dtype=np.float64)
    obs_total_bases = np.nansum(dp_corrected * bin_sizes[:, None], axis=0)
    bb_lengths = np.bincount(
        bin_df["bb_id"].to_numpy(np.int64), weights=bin_sizes, minlength=num_bbs
    )

    for o in tumor_dataset_indices:
        m = get_rdr_base_dataset_id.get(o)
        if m is not None:
            lib = obs_total_bases[m] / obs_total_bases[o]
            logging.info(
                f"  bb RDR {dataset_ids[o]} / base {dataset_ids[m]}: library factor={lib:.4f}"
            )
            base = bb_dp[:, m]
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[:, rdr_pos[o]] = np.divide(
                    bb_dp[:, o] * lib,
                    base,
                    where=base > 0,
                    out=np.full(num_bbs, np.nan, dtype=np.float32),
                )
        else:
            vals = bb_dp[:, o]
            positive = np.isfinite(vals) & (vals > 0)
            if positive.any():
                med = np.median(vals[positive])
                logging.info(
                    f"  bb median-centering {dataset_ids[o]}: median={med:.4f}"
                )
                finite = np.isfinite(vals)
                bb_rdr[finite, rdr_pos[o]] = vals[finite] / med
            else:
                logging.warning(
                    f"  bb median-centering {dataset_ids[o]}: no positive finite bb, "
                    "the whole RDR column stays NaN"
                )
        nan_col = np.isnan(bb_rdr[:, rdr_pos[o]])
        log_ratios(
            f"bb RDR NaN {dataset_ids[o]}",
            int(nan_col.sum()),
            num_bbs,
            bb_lengths[nan_col],
            "bb",
            prefix="  ",
        )

    return bb_rdr
