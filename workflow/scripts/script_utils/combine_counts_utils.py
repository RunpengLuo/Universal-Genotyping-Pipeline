"""Allele/SNP helpers shared by phase_and_concat and combine_counts.

Five groups, in pipeline order:

1. allele matrices - per-replicate counts onto the parent SNP set, then stacked
2. observations    - which matrix column belongs to which replicate/assay
3. SNP set         - the union SNP table and the per-SNP range split
4. SNP filters     - one bool mask each, combined by ``apply_masks_to_df``
5. depth and RDR   - fixed-bin depth onto bbs, then the RDR ratio
"""

import logging

import numpy as np
import pandas as pd

from scipy.sparse import csr_matrix, hstack
from scipy.stats import beta

from io_utils import read_BED
from range_utils import assign_pos_to_range, assign_range_to_range, overlaps_any_range
from utils import sort_df_chr


##################################################
# allele matrices: per-replicate counts -> one joint matrix


def map_allele_mat_to_snps(
    parent_keys: pd.Index,
    child_snps: pd.DataFrame,
    tot_mtx: csr_matrix,
    ad_mtx: csr_matrix,
    ncells: int,
):
    """Map a replicate's DP/AD matrices onto the parent SNP set (0-fill missing loci).

    Args:
        parent_keys: parent SNP KEYs (``#CHROM_POS``) in output row order.
        child_snps: DataFrame with ``KEY`` + ``RAW_SNP_DF_IDX`` for this replicate's loci
            (in the matrices' row order).
        tot_mtx, ad_mtx: ``(m, ncells)`` count matrices, ``m == len(child_snps)``.
        ncells: number of observations (1 for bulk pseudobulk; #barcodes for single-cell).

    Returns:
        ``(tot_canon, ad_canon)`` csr matrices of shape ``(len(parent_keys), ncells)``.
    """
    M = len(parent_keys)
    m = len(child_snps)

    tot_mtx = tot_mtx.tocsr()
    ad_mtx = ad_mtx.tocsr()
    assert tot_mtx.shape == ad_mtx.shape, (
        f"count matrices, shape mismatch: {tot_mtx.shape} vs {ad_mtx.shape}"
    )
    assert tot_mtx.shape == (m, ncells), (
        f"count matrix, shape {tot_mtx.shape}, expected {(m, ncells)}"
    )

    raw_snp_df_idx = child_snps["RAW_SNP_DF_IDX"].to_numpy()
    tot_mtx = tot_mtx[raw_snp_df_idx, :]
    ad_mtx = ad_mtx[raw_snp_df_idx, :]

    dup_mask = child_snps["KEY"].duplicated(keep=False).to_numpy()
    n_dup_snps = int(dup_mask.sum())
    if n_dup_snps > 0:
        n_dup_keys = int(child_snps.loc[dup_mask, "KEY"].nunique())
        logging.warning(
            f"#found duplicated SNP #CHR/POS, {n_dup_snps}/{m} SNPs over {n_dup_keys} keys, drop all"
        )
        child_snps = child_snps.loc[~dup_mask].reset_index(drop=True)
        tot_mtx = tot_mtx[~dup_mask, :]
        ad_mtx = ad_mtx[~dup_mask, :]

    m = len(child_snps)
    child_keys = pd.Index(child_snps["KEY"])

    child_loc = child_keys.get_indexer(parent_keys)
    present = child_loc >= 0
    logging.info(f"located SNPs in parent={present.sum()}/{M}")

    common_pidx = np.flatnonzero(present)
    common_cidx = child_loc[present]

    tot_present = tot_mtx[common_cidx, :].tocoo()
    tot_canon = csr_matrix(
        (tot_present.data, (common_pidx[tot_present.row], tot_present.col)),
        shape=(M, ncells),
        dtype=tot_mtx.dtype,
    )

    ad_present = ad_mtx[common_cidx, :].tocoo()
    ad_canon = csr_matrix(
        (ad_present.data, (common_pidx[ad_present.row], ad_present.col)),
        shape=(M, ncells),
        dtype=ad_mtx.dtype,
    )

    return tot_canon, ad_canon


def hstack_replicate_mats(tot_list: list, ad_list: list):
    """Horizontally stack per-replicate total and alt-count matrices.

    Derives the REF matrix as ``TOT - ALT``.

    Parameters
    ----------
    tot_list : list of csr_matrix
        Per-replicate total depth matrices (SNP features x observations).
    ad_list : list of csr_matrix
        Per-replicate alt-allele count matrices (SNP features x observations).

    Returns
    -------
    tuple[csr_matrix, csr_matrix, csr_matrix]
        ``(tot_mtx, ref_mtx, alt_mtx)`` concatenated across replicates.
    """
    if len(tot_list) == 1:
        return tot_list[0], tot_list[0] - ad_list[0], ad_list[0]
    tot_mtx = hstack(tot_list, format="csr")
    alt_mtx = hstack(ad_list, format="csr")
    ref_mtx = tot_mtx - alt_mtx
    return tot_mtx, ref_mtx, alt_mtx


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
# SNP set: the union table and the per-SNP range split


def build_union_snps(snps_list):
    """Union per-assay SNP tables into one genomically-sorted set.

    Keeps the shared annotation columns (``PS``/``feature_id`` only when present in
    EVERY assay), dedupes on ``(#CHR, POS0)``, sorts, and adds a 0-based
    ``snp_id``. Returns ``(snps, has_ps, has_feature)``.
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
    return snps, has_ps, has_feature


def build_pos_ranges(snps: pd.DataFrame, regions: pd.DataFrame, colname="region_id"):
    """Split each region into one ``[START, END)`` range per position it contains.

    A position owns the span between the midpoints of its neighbours, bounded by the
    region edges; a lone position owns the whole region. Positions outside every region
    keep ``START = END = 0`` and an empty ``region_id``. Membership is by 0-based
    ``POS0``, like every other range operation.

    Args:
        snps: Positions with ``#CHR``, ``POS0``, sorted genomically.
        regions: Ranges with ``#CHR``, ``START``, ``END``, ``region_id`` and
            optionally ``seg_id``.
        colname: Column to write the region identifier into.

    Returns:
        *snps* with ``START``, ``END``, ``BLOCKSIZE``, *colname* (and ``seg_id``).
    """
    has_seg = "seg_id" in regions.columns
    regions = regions.reset_index(drop=True)
    regions["_reg_id"] = np.arange(len(regions))

    snps, _ = assign_pos_to_range(snps, regions, ref_id="_reg_id", pos_col="POS0")
    snps["START"] = 0
    snps["END"] = 0
    snps[colname] = ""
    if has_seg:
        snps["seg_id"] = ""

    inside = snps["_reg_id"].notna()
    n_out = int((~inside).sum())
    if n_out:
        logging.info(f"SNP ranges: {n_out}/{len(snps)} SNPs outside every region")

    reg_start = regions["START"].to_numpy()
    reg_end = regions["END"].to_numpy()
    for reg_id, grp in snps.loc[inside].groupby("_reg_id", sort=False):
        idx = grp.index.to_numpy()
        r = int(reg_id)
        snps.loc[idx, colname] = regions.at[r, "region_id"]
        if has_seg:
            snps.loc[idx, "seg_id"] = regions.at[r, "seg_id"]
        pos0 = grp["POS0"].to_numpy()
        if len(idx) == 1:
            bounds = np.array([reg_start[r], reg_end[r]])
        else:
            mids = np.ceil((pos0[:-1] + pos0[1:]) / 2).astype(np.uint32)
            bounds = np.concatenate([[reg_start[r]], mids, [reg_end[r]]])
        snps.loc[idx, "START"] = bounds[:-1]
        snps.loc[idx, "END"] = bounds[1:]

    snps.drop(columns="_reg_id", inplace=True)
    snps["BLOCKSIZE"] = snps["END"] - snps["START"]
    return snps


##################################################
# SNP filters: one bool mask each, combined by apply_masks_to_df


def apply_masks_to_df(df: pd.DataFrame, *masks):
    """Keep the rows of *df* that every mask keeps, reindexed from 0.

    Args:
        df: Frame to filter. Copied, not modified.
        *masks: One or more bool arrays of length ``len(df)``.

    Returns:
        ``(df, keep)`` - the surviving rows, ``reset_index(drop=True)``, and the
        combined mask over the INPUT rows, for subsetting a parallel matrix the
        same way.
    """
    assert masks, "apply_masks_to_df, no mask given"
    keep = np.logical_and.reduce([np.asarray(m, dtype=bool) for m in masks])
    assert len(keep) == len(df), f"mask length {len(keep)} != {len(df)} rows"
    return df.loc[keep].reset_index(drop=True), keep


def get_mask_by_region(snps: pd.DataFrame, regions: pd.DataFrame):
    """Keep SNPs inside any region."""
    mask = overlaps_any_range(snps, regions)
    logging.info(f"filter by region, #passed SNPs={np.sum(mask)}/{len(snps)}")
    return mask


def get_mask_by_blacklist(snps: pd.DataFrame, blacklist_bed):
    """Drop SNPs inside any blacklisted range; all-True when no blacklist is given."""
    if blacklist_bed is None:
        return np.ones(len(snps), dtype=bool)
    hit = overlaps_any_range(snps, read_BED(blacklist_bed))
    num_passes = len(snps) - np.sum(hit)
    logging.info(f"filter by blacklist, #passed SNPs={num_passes}/{len(snps)}")
    return ~hit


def get_mask_by_exon(snps: pd.DataFrame):
    """Keep only exonic SNPs.

    Reads the ``feature_type`` column that ``feature_utils.annotate_feature_type`` writes.
    The caller gates the filter on the ``exon_only`` config key.
    """
    is_exon = (snps["feature_type"] == "exon").to_numpy()
    logging.info(f"filter by exon, #passed SNPs={np.sum(is_exon)}/{len(snps)}")
    return is_exon


def get_mask_by_depth(snps: pd.DataFrame, tot_mtx: csr_matrix, min_dp=1):
    """Return a boolean mask keeping SNPs where every sample meets the minimum depth.

    Parameters
    ----------
    snps : pd.DataFrame
        SNP info DataFrame (used only for logging).
    tot_mtx : csr_matrix
        Total depth matrix (SNPs x samples).
    min_dp : int
        Minimum depth threshold per sample.

    Returns
    -------
    np.ndarray
        Boolean mask of length ``len(snps)``.
    """
    mask = np.all(tot_mtx >= min_dp, axis=1)
    logging.info(
        f"filter by depth, min_dp={min_dp}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


def get_mask_by_het_balanced(
    snps: pd.DataFrame,
    ref_mtx: csr_matrix,
    alt_mtx: csr_matrix,
    gamma: float,
    normal_idx=0,
):
    """
    mask SNPs if normal sample failed beta-posterior credible interval test with beta(1, 1) prior.
    """
    p_lower = gamma / 2.0
    p_upper = 1.0 - p_lower
    q = np.array([p_lower, p_upper])
    het_cred_ints = beta.ppf(
        q[None, :],
        ref_mtx[:, normal_idx][:, None] + 1,
        alt_mtx[:, normal_idx][:, None] + 1,
    )
    mask = (het_cred_ints[:, 0] <= 0.5) & (0.5 <= het_cred_ints[:, 1])
    logging.info(
        f"filter by balanced Het-SNPs on normal sample, gamma={gamma}, #passed SNPs={np.sum(mask)}/{len(snps)}"
    )
    return mask


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
