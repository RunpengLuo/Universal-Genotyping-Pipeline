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

from range_utils import assign_range_to_range


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


def summarize_read_depth_bb(
    assay2dataset_indices, bin_df, dp_bin_dfs, dp_corrected_list, num_bbs, num_datasets
):
    """Length-weighted aggregation of corrected fixed-bin depth into bbs.

    Two bin frames, on purpose: *bin_df* is the union grid carrying ``bb_id``, and each
    *dp_bin_dfs* entry is one assay's own fixed bins, row-aligned to its
    *dp_corrected_list* matrix. Each assay's bins are assigned to a bb by midpoint over
    *bin_df*, so a finer WES bin set projects onto the WGS bbs.

    Args:
        assay2dataset_indices: ``{assay_type: dataset indices}``, in the order of
            *dp_bin_dfs* and *dp_corrected_list*. Column ``s`` of an assay's matrix goes
            to its ``s``-th index, so one assay's datasets need not be adjacent.

    Returns:
        ``(bb_dp, bb_bases)``: per-bb mean depth and per-bb total aligned bases, both
        ``(num_bbs, num_datasets)``.
    """

    bb_spans = (
        bin_df.groupby("bb_id", sort=True)
        .agg(
            **{
                "#CHR": ("#CHR", "first"),
                "START": ("START", "min"),
                "END": ("END", "max"),
            }
        )
        .reset_index()
    )
    bb_dp = np.full((num_bbs, num_datasets), np.nan, dtype=np.float32)
    bb_bases = np.zeros((num_bbs, num_datasets), dtype=np.float64)
    for (at, dataset_indices), bins_a, dp_a in zip(
        assay2dataset_indices.items(), dp_bin_dfs, dp_corrected_list
    ):
        # a finer WES bin set projects onto the WGS bbs by midpoint
        mapped, na_idx = assign_range_to_range(
            bins_a[["#CHR", "START", "END"]], bb_spans, "bb_id", rule="midpoint"
        )
        bb_ids = mapped["bb_id"].fillna(-1).to_numpy(np.int64)
        valid = bb_ids >= 0
        if len(na_idx):
            logging.info(
                f"{at}: {len(na_idx)}/{len(bins_a)} fixed bins outside all bbs (dropped)"
            )
        vb = bb_ids[valid]
        bin_lengths = (bins_a["END"] - bins_a["START"]).to_numpy(dtype=np.float64)[
            valid
        ]
        total_len_per_bb = np.bincount(vb, weights=bin_lengths, minlength=num_bbs)
        for s, obs in enumerate(dataset_indices):
            weighted_sums = np.bincount(
                vb, weights=dp_a[valid, s] * bin_lengths, minlength=num_bbs
            )
            bb_bases[:, obs] = weighted_sums
            with np.errstate(invalid="ignore"):
                bb_dp[:, obs] = weighted_sums / total_len_per_bb
    return bb_dp, bb_bases


def summarize_rdr_bb(
    assay2dataset_indices,
    dp_bin_dfs,
    dp_corrected_list,
    bb_dp,
    tumor_dataset_indices,
    get_rdr_base_dataset_id,
    rdr_outlier_quantile,
    dataset_ids,
):
    """Per-bb RDR for every tumor observation.

    Each tumor with an RDR base observation (``{tumor index: base index}``) is normalized
    by that base, library-size corrected; a tumor without a base is median-centered. The base
    may be any observation (e.g. a different assay/platform), so library sizes are
    computed globally per observation. Entries above the ``1 - rdr_outlier_quantile``
    quantile are set to NaN. Returns a ``(num_bbs, len(tumor_dataset_indices))`` array
    aligned to *tumor_dataset_indices*; *assay2dataset_indices* is as in
    ``summarize_read_depth_bb``.
    """
    num_bbs, num_datasets = bb_dp.shape
    bb_rdr = np.full((num_bbs, len(tumor_dataset_indices)), np.nan, dtype=np.float32)
    rdr_pos = {o: i for i, o in enumerate(tumor_dataset_indices)}

    # global per-observation total aligned bases for library-size correction
    obs_total_bases = np.full(num_datasets, np.nan, dtype=np.float64)
    for dataset_indices, bins_a, dp_a in zip(
        assay2dataset_indices.values(), dp_bin_dfs, dp_corrected_list
    ):
        bin_sizes = (bins_a["END"] - bins_a["START"]).to_numpy(dtype=np.float64)
        obs_total_bases[dataset_indices] = np.nansum(dp_a * bin_sizes[:, None], axis=0)

    for o in tumor_dataset_indices:
        m = get_rdr_base_dataset_id.get(o)
        if m is not None:
            lib = obs_total_bases[m] / obs_total_bases[o]
            logging.info(
                f"  bb RDR {dataset_ids[o]} / base {dataset_ids[m]}: library factor={lib:.4f}"
            )
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[:, rdr_pos[o]] = bb_dp[:, o] / bb_dp[:, m] * lib
        else:
            vals = bb_dp[:, o]
            valid_i = np.isfinite(vals) & (vals > 0)
            if valid_i.any():
                med = np.median(vals[valid_i])
                logging.info(
                    f"  bb median-centering {dataset_ids[o]}: median={med:.4f}"
                )
                with np.errstate(invalid="ignore", divide="ignore"):
                    bb_rdr[valid_i, rdr_pos[o]] = vals[valid_i] / med

    if rdr_outlier_quantile > 0:
        rdr_upper = np.nanquantile(bb_rdr, 1 - rdr_outlier_quantile)
        n_outlier = int(np.nansum(bb_rdr > rdr_upper))
        logging.info(
            f"RDR outlier filter: quantile={rdr_outlier_quantile}, "
            f"threshold={rdr_upper:.4f}, {n_outlier} entries set to NaN"
        )
        bb_rdr[bb_rdr > rdr_upper] = np.nan
    return bb_rdr
