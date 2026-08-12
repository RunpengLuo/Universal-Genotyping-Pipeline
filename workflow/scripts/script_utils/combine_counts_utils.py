"""Observation bookkeeping and bulk depth/RDR summaries for combine_counts.

Last update: 2026-08-11

Functions:
- observation_cluster_ids: map each matrix column to its roster row
- summarize_read_depth_bb: length-weighted fixed-bin depth aggregated per bb
- summarize_rdr_bb: per-bb RDR, matched-normal or median-centered
"""

import logging

import numpy as np
import pandas as pd


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


##################################################
# bulk read depth and RDR


def summarize_read_depth_bb(bin_df, dp_corrected, num_bbs):
    """Length-weighted aggregation of corrected fixed-bin depth into bbs.

    *bin_df* is the shared fixed-bin grid carrying the ``bb_id`` that
    ``build_adaptive_bins`` stamped, row-aligned to *dp_corrected*. NaN is masked per
    column: a window the correction could not resolve for one dataset still counts for
    every other dataset.

    Args:
        bin_df: Fixed bins with ``START``, ``END`` and ``bb_id``.
        dp_corrected: ``(n_bins, n_datasets)`` corrected depth, NaN where undefined.
        num_bbs: Number of bbs.

    Returns:
        ``(bb_dp, bb_bases)``: per-bb mean depth and per-bb total aligned bases, both
        ``(num_bbs, n_datasets)``; *bb_dp* is NaN for a bb with no finite window.
    """
    bb_ids = bin_df["bb_id"].to_numpy(np.int64)
    bin_lengths = (bin_df["END"] - bin_df["START"]).to_numpy(dtype=np.float64)
    num_datasets = dp_corrected.shape[1]
    bb_dp = np.full((num_bbs, num_datasets), np.nan, dtype=np.float32)
    bb_bases = np.zeros((num_bbs, num_datasets), dtype=np.float64)
    for s in range(num_datasets):
        finite = np.isfinite(dp_corrected[:, s])
        n_nan = len(finite) - int(finite.sum())
        if n_nan:
            logging.info(f"  dataset column {s}: {n_nan} NaN fixed bins skipped")
        ids, lengths = bb_ids[finite], bin_lengths[finite]
        weighted_sums = np.bincount(
            ids, weights=dp_corrected[finite, s] * lengths, minlength=num_bbs
        )
        total_len_per_bb = np.bincount(ids, weights=lengths, minlength=num_bbs)
        bb_bases[:, s] = weighted_sums
        with np.errstate(invalid="ignore", divide="ignore"):
            bb_dp[:, s] = weighted_sums / total_len_per_bb
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

    return bb_rdr
