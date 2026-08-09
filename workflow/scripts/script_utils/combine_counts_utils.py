"""Observation, SNP-set and RDR helpers for the combine_counts scripts.

Three groups, in pipeline order:

1. observations  - which matrix column belongs to which replicate/assay
2. SNP set       - the union SNP table across assays
3. depth and RDR - fixed-bin depth onto bbs, then the RDR ratio

The upstream half of the pipeline is ``phase_and_concat_utils``; nothing is shared.
"""

import logging

import numpy as np
import pandas as pd

from range_utils import assign_range_to_range
from utils import sort_df_chr


##################################################
# observations: which matrix column is which replicate/assay


def observation_cluster_ids(rep2bc: pd.DataFrame, dataset_ids):
    """Cluster id per observation: REP_ID,BARCODE -> int64 index into dataset_ids.

    Categorical mapping with explicit ``dataset_ids`` order ensures the codes
    align with the position of each dataset_id in the caller's dataset_ids list.
    """
    cats = pd.Categorical(rep2bc["REP_ID"], categories=list(dataset_ids))
    codes = np.asarray(cats.codes, dtype=np.int64)
    assert (codes >= 0).all(), "barcodes.full, REP_ID values outside dataset_ids"
    return codes


def build_assay_obs_clusters(sample_df, bulk_assays):
    """Cluster the joint observations by assay into per-assay descriptors.

    Each cluster records the assay's observation ``offset``, size ``n``, and
    ``tumor_obs``. Returns ``(assay_obs_clusters, tumor_obs_all)``.
    """
    dataset_assays = sample_df["assay_type"].tolist()
    sample_types = sample_df["sample_type"].tolist()
    assay_obs_clusters = []
    for at in bulk_assays:
        obs = [i for i, a in enumerate(dataset_assays) if a == at]
        assert obs, f"joint sample sheet, no sample for assay {at}"
        stypes = [sample_types[i] for i in obs]
        assay_obs_clusters.append(
            {
                "assay": at,
                "offset": obs[0],
                "n": len(obs),
                "tumor_obs": [obs[i] for i, st in enumerate(stypes) if st == "tumor"],
            }
        )
    tumor_obs_all = [o for c in assay_obs_clusters for o in c["tumor_obs"]]
    return assay_obs_clusters, tumor_obs_all


##################################################
# SNP set: the union table across assays


def build_union_snps(snps_list):
    """Union per-assay SNP tables into one genomically-sorted set.

    Keeps the shared annotation columns (``PS``/``feature_id`` only when present in
    EVERY assay), dedupes on ``(#CHR, POS0)``, sorts, and adds a 0-based
    ``snp_id``. Returns ``(snps, has_feature)``; ``setup_phase_clusters`` reports on
    ``PS``.
    """
    has_ps = all("PS" in s.columns for s in snps_list)
    has_feature = all("feature_id" in s.columns for s in snps_list)
    has_seg = all("seg_id" in s.columns for s in snps_list)
    annot_cols = (
        ["#CHR", "POS", "POS0", "START", "END", "region_id"]
        + (["seg_id"] if has_seg else [])
        + (["PS"] if has_ps else [])
        + (["feature_id"] if has_feature else [])
    )
    # feature_id is GTF-derived per assay, so the same SNP carries an identical
    # string in every assay -> keep-first dedup is already a correct cross-assay union
    snps = pd.concat(
        [s[annot_cols] for s in snps_list], ignore_index=True
    ).drop_duplicates(["#CHR", "POS0"])
    snps = sort_df_chr(snps, ch="#CHR", pos="POS0").reset_index(drop=True)
    snps["snp_id"] = np.arange(len(snps))
    return snps, has_feature


##################################################
# bulk read depth and RDR


def aggregate_bin_depth_to_bbs(
    assay_obs_clusters, scaffold, bin_df_list, dp_corrected_list, num_bbs, total_samples
):
    """Length-weighted aggregation of corrected fixed-bin depth into bbs.

    Each assay's fixed bins are assigned to a bb by midpoint over the ``scaffold``
    (the binning fixed bins carrying ``bb_id``), so a finer WES bin set is projected onto
    the WGS bbs. Returns ``(bb_dp, bb_bases)``: per-bb mean depth and per-bb total
    aligned bases, both ``(num_bbs, total_samples)``.
    """

    bb_spans = (
        scaffold.groupby("bb_id", sort=True)
        .agg(
            **{
                "#CHR": ("#CHR", "first"),
                "START": ("START", "min"),
                "END": ("END", "max"),
            }
        )
        .reset_index()
    )
    bb_dp = np.full((num_bbs, total_samples), np.nan, dtype=np.float32)
    bb_bases = np.zeros((num_bbs, total_samples), dtype=np.float64)
    for clu, bins_a, dp_a in zip(assay_obs_clusters, bin_df_list, dp_corrected_list):
        # a finer WES bin set projects onto the WGS bbs by midpoint
        mapped, na_idx = assign_range_to_range(
            bins_a[["#CHR", "START", "END"]], bb_spans, "bb_id", rule="midpoint"
        )
        bb_ids = mapped["bb_id"].fillna(-1).to_numpy(np.int64)
        valid = bb_ids >= 0
        if len(na_idx):
            logging.info(
                f"{clu['assay']}: {len(na_idx)}/{len(bins_a)} fixed bins outside all bbs (dropped)"
            )
        vb = bb_ids[valid]
        bin_lengths = (bins_a["END"] - bins_a["START"]).to_numpy(dtype=np.float64)[
            valid
        ]
        total_len_per_bb = np.bincount(vb, weights=bin_lengths, minlength=num_bbs)
        for s in range(clu["n"]):
            weighted_sums = np.bincount(
                vb, weights=dp_a[valid, s] * bin_lengths, minlength=num_bbs
            )
            bb_bases[:, clu["offset"] + s] = weighted_sums
            with np.errstate(invalid="ignore"):
                bb_dp[:, clu["offset"] + s] = weighted_sums / total_len_per_bb
    return bb_dp, bb_bases


def build_rdr_base_map(sample_df):
    """Map each tumor observation to its RDR base (denominator) observation.

    A tumor row's optional ``RDR_BASE_REP_ID`` names the ``REP_ID`` of the
    sample used as its RDR baseline. Returns ``{tumor_obs: base_obs}``; a tumor
    with an unset ``RDR_BASE_REP_ID`` is omitted (median-normalized downstream).
    """
    dataset_ids = sample_df["REP_ID"].tolist()
    sample_types = sample_df["sample_type"].tolist()
    dataset_id_to_obs = {rid: i for i, rid in enumerate(dataset_ids)}

    has_col = "RDR_BASE_REP_ID" in sample_df.columns
    base_dataset_ids = (
        sample_df["RDR_BASE_REP_ID"].tolist() if has_col else [None] * len(dataset_ids)
    )

    base_map = {}
    for i in range(len(dataset_ids)):
        if sample_types[i] != "tumor":
            continue
        base_dataset_id = base_dataset_ids[i]
        if has_col and pd.notna(base_dataset_id) and str(base_dataset_id) != "":
            assert base_dataset_id in dataset_id_to_obs, (
                f"{dataset_ids[i]}: RDR_BASE_REP_ID {base_dataset_id!r} is not a REP_ID"
            )
            assert dataset_id_to_obs[base_dataset_id] != i, (
                f"{dataset_ids[i]}: RDR_BASE_REP_ID {base_dataset_id!r} is itself"
            )
            base_map[i] = dataset_id_to_obs[base_dataset_id]
    return base_map


def compute_bb_rdr(
    assay_obs_clusters,
    bin_df_list,
    dp_corrected_list,
    bb_dp,
    tumor_obs_all,
    base_map,
    rdr_outlier_quantile,
    dataset_ids,
):
    """Per-bb RDR for every tumor observation.

    Each tumor with an RDR base observation (from ``base_map``) is normalized by that
    base, library-size corrected; a tumor without a base is median-centered. The base
    may be any observation (e.g. a different assay/platform), so library sizes are
    computed globally per observation. Entries above the ``1 - rdr_outlier_quantile``
    quantile are set to NaN. Returns a ``(num_bbs, len(tumor_obs_all))`` array aligned
    to ``tumor_obs_all``.
    """
    num_bbs, total_samples = bb_dp.shape
    bb_rdr = np.full((num_bbs, len(tumor_obs_all)), np.nan, dtype=np.float32)
    rdr_pos = {o: i for i, o in enumerate(tumor_obs_all)}

    # global per-observation total aligned bases for library-size correction
    obs_total_bases = np.full(total_samples, np.nan, dtype=np.float64)
    for clu, bins_a, dp_a in zip(assay_obs_clusters, bin_df_list, dp_corrected_list):
        bin_sizes = (bins_a["END"] - bins_a["START"]).to_numpy(dtype=np.float64)
        tb = np.nansum(dp_a * bin_sizes[:, None], axis=0)
        for s in range(clu["n"]):
            obs_total_bases[clu["offset"] + s] = tb[s]

    for o in tumor_obs_all:
        m = base_map.get(o)
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
