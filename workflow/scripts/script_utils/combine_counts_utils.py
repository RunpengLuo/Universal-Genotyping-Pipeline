"""Utility functions for SNP/allele processing, phasing, and filtering.

Used by phase_and_concat_{bulk,single_cell}.py, combine_counts.py, and combine_counts_nonbulk.py.
"""

import logging

import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy.sparse import csr_matrix, hstack
from scipy.stats import beta

from interval_utils import assign_pos_to_range
from io_utils import read_VCF
from utils import sort_df_chr


def canon_mat_from_files(
    parent_keys: pd.Index,
    vcf_file: str,
    tot_mtx_file: str,
    ad_mtx_file: str,
    ncells: int,
):
    """Read a replicate's cellSNP-format files, then canonicalize onto the parent grid.

    Thin file-reading wrapper over ``canon_mat_one_replicate`` for the cellsnp-lite path
    (single-cell). Bulk builds the DataFrame/matrices from bcftools counts instead.
    """
    child_snps = read_VCF(vcf_file, addkey=True)
    tot_mtx = mmread(tot_mtx_file).tocsr()
    ad_mtx = mmread(ad_mtx_file).tocsr()
    return canon_mat_one_replicate(parent_keys, child_snps, tot_mtx, ad_mtx, ncells)


def bcftools_counts_to_child_mats(bcf_df: pd.DataFrame, parent_alt_by_key: dict):
    """Convert a bcftools counts DataFrame into ``(child_snps, tot_mtx, ad_mtx)`` for canon.

    ``ALT`` per locus is the depth of the *parent* ALT allele (matched against the bcftools
    ALT list; 0 when the parent ALT was not observed), ``DP = ref + alt`` (usable het depth,
    so downstream ``REF = DP - ALT`` is exact). Rows follow ``bcf_df`` order.

    Args:
        bcf_df: output of ``read_bcftools_counts`` (KEY, ALT list, AD list, RAW_SNP_IDX).
        parent_alt_by_key: parent ALT allele keyed by ``#CHROM_POS``.

    Returns:
        ``(child_snps, tot_mtx, ad_mtx)``: child_snps has KEY + RAW_SNP_IDX; matrices are
        ``(len(bcf_df), 1)`` csr.
    """
    n = len(bcf_df)
    keys = bcf_df["KEY"].to_numpy()
    alt_lists = bcf_df["ALT"].to_numpy()
    ad_lists = bcf_df["AD"].to_numpy()
    tot = np.zeros(n, dtype=np.int64)
    ad = np.zeros(n, dtype=np.int64)
    for i in range(n):
        adv = ad_lists[i]
        ref_cnt = adv[0] if len(adv) else 0
        palt = parent_alt_by_key.get(keys[i])
        alts = alt_lists[i]
        alt_cnt = (
            adv[1 + alts.index(palt)] if (palt is not None and palt in alts) else 0
        )
        ad[i] = alt_cnt
        tot[i] = ref_cnt + alt_cnt
    child_snps = pd.DataFrame({"KEY": keys, "RAW_SNP_IDX": np.arange(n)})
    tr = np.flatnonzero(tot)
    ar = np.flatnonzero(ad)
    tot_mtx = csr_matrix(
        (tot[tr], (tr, np.zeros(len(tr), dtype=np.int64))), shape=(n, 1)
    )
    ad_mtx = csr_matrix((ad[ar], (ar, np.zeros(len(ar), dtype=np.int64))), shape=(n, 1))
    return child_snps, tot_mtx, ad_mtx


def canon_mat_one_replicate(
    parent_keys: pd.Index,
    child_snps: pd.DataFrame,
    tot_mtx: csr_matrix,
    ad_mtx: csr_matrix,
    ncells: int,
):
    """Map a replicate's DP/AD matrices onto the parent SNP grid (0-fill missing loci).

    Args:
        parent_keys: parent SNP KEYs (``#CHROM_POS``) in output row order.
        child_snps: DataFrame with ``KEY`` + ``RAW_SNP_IDX`` for this replicate's loci
            (in the matrices' row order).
        tot_mtx, ad_mtx: ``(m, ncells)`` count matrices, ``m == len(child_snps)``.
        ncells: number of columns (1 for bulk pseudobulk; #barcodes for single-cell).

    Returns:
        ``(tot_canon, ad_canon)`` csr matrices of shape ``(len(parent_keys), ncells)``.
    """
    M = len(parent_keys)
    m = len(child_snps)

    tot_mtx = tot_mtx.tocsr()
    ad_mtx = ad_mtx.tocsr()
    assert tot_mtx.shape == ad_mtx.shape
    assert tot_mtx.shape == (m, ncells)

    raw_snp_ids = child_snps["RAW_SNP_IDX"].to_numpy()
    tot_mtx = tot_mtx[raw_snp_ids, :]
    ad_mtx = ad_mtx[raw_snp_ids, :]

    dup_mask = child_snps["KEY"].duplicated(keep=False).to_numpy()
    n_dup_rows = int(dup_mask.sum())
    if n_dup_rows > 0:
        n_dup_keys = int(child_snps.loc[dup_mask, "KEY"].nunique())
        logging.warning(
            f"#found duplicated SNP #CHR/POS, #rows={n_dup_rows}/{m} #SNPs={n_dup_keys}, drop all"
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


def cell_rep_idx_from_mapping(rep2bc: pd.DataFrame, dataset_ids):
    """Convert a REP_ID,BARCODE DataFrame into an int64 array of rep indices.

    Categorical mapping with explicit ``dataset_ids`` order ensures the codes
    align with the position of each rep in the caller's dataset_ids list.
    """
    cats = pd.Categorical(rep2bc["REP_ID"], categories=list(dataset_ids))
    codes = np.asarray(cats.codes, dtype=np.int64)
    assert (codes >= 0).all(), (
        "barcodes.full contains REP_ID values outside dataset_ids"
    )
    return codes


def merge_mats(tot_list: list, ad_list: list):
    """Horizontally stack per-replicate total and alt-count matrices.

    Derives the REF matrix as ``TOT - ALT``.

    Parameters
    ----------
    tot_list : list of csr_matrix
        Per-replicate total depth matrices (SNPs x cells).
    ad_list : list of csr_matrix
        Per-replicate alt-allele count matrices (SNPs x cells).

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
# combine_counts: SNP parsing, grouping, and bulk depth/RDR aggregation


def build_union_snp_grid(snps_list):
    """Union per-assay SNP tables onto one genomically-sorted grid.

    Keeps shared annotation columns (``PS``/``feature_id`` only when present in
    EVERY assay), dedupes on ``(#CHR, POS0)``, sorts, and adds a 0-based
    ``snp_row``. Returns ``(snps, has_ps, has_feature)``.
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
    snps["snp_row"] = np.arange(len(snps))
    return snps, has_ps, has_feature


def build_assay_blocks(sample_df, bulk_assays):
    """Group joint sample columns by assay into per-assay block descriptors.

    Each block records the assay's column ``offset``, size ``n``, and
    ``tumor_cols``. Returns ``(assay_blocks, tumor_cols_all)``.
    """
    col_assay = sample_df["assay_type"].tolist()
    col_stype = sample_df["sample_type"].tolist()
    assay_blocks = []
    for at in bulk_assays:
        cols = [i for i, a in enumerate(col_assay) if a == at]
        assert cols, f"no samples for assay {at} in joint sample sheet"
        stypes = [col_stype[i] for i in cols]
        assay_blocks.append(
            {
                "assay": at,
                "offset": cols[0],
                "n": len(cols),
                "tumor_cols": [cols[i] for i, st in enumerate(stypes) if st == "tumor"],
            }
        )
    tumor_cols_all = [c for blk in assay_blocks for c in blk["tumor_cols"]]
    return assay_blocks, tumor_cols_all


def build_rdr_base_map(sample_df):
    """Map each tumor column to its RDR base (denominator) column.

    A tumor row's optional ``RDR_BASE_REP_ID`` names the ``REP_ID`` of the
    sample used as its RDR baseline. Returns ``{tumor_col: base_col}``; a tumor
    with an unset ``RDR_BASE_REP_ID`` is omitted (median-normalized downstream).
    """
    col_repid = sample_df["REP_ID"].tolist()
    col_stype = sample_df["sample_type"].tolist()
    repid_to_col = {rid: i for i, rid in enumerate(col_repid)}

    has_col = "RDR_BASE_REP_ID" in sample_df.columns
    col_base_rep = (
        sample_df["RDR_BASE_REP_ID"].tolist() if has_col else [None] * len(col_repid)
    )

    base_map = {}
    for i in range(len(col_repid)):
        if col_stype[i] != "tumor":
            continue
        brep = col_base_rep[i]
        if has_col and pd.notna(brep) and str(brep) != "":
            assert brep in repid_to_col, (
                f"RDR_BASE_REP_ID={brep!r} for tumor {col_repid[i]!r} is not a REP_ID"
            )
            assert repid_to_col[brep] != i, (
                f"RDR_BASE_REP_ID={brep!r} for tumor {col_repid[i]!r} is itself"
            )
            base_map[i] = repid_to_col[brep]
    return base_map


def aggregate_window_depth_to_bins(
    assay_blocks, scaffold, window_df_list, dp_corrected_list, num_bbs, total_samples
):
    """Length-weighted aggregation of corrected window depth into adaptive bins.

    Each assay's windows are assigned to bins by midpoint over the ``scaffold``
    (the binning-grid windows carrying ``bin_id``), so a finer WES grid is
    projected onto the WGS bins. Returns ``(bb_dp, bb_bases)``: per-bin mean depth
    and per-bin total aligned bases, both ``(num_bbs, total_samples)``.
    """

    def _windows_to_bins(win_a, bin_spans):
        """Assign each window to the bin whose genomic span contains its midpoint.

        Returns an int64 array of ``bin_id`` per window, ``-1`` where the midpoint
        falls in no bin. Works for a same-grid assay (window in its own bin) and for
        a finer WES grid projected onto WGS bins.
        """
        mids = pd.DataFrame(
            {
                "#CHR": win_a["#CHR"].to_numpy(),
                "POS0": ((win_a["START"] + win_a["END"]) // 2).to_numpy(np.int64),
            }
        )
        mids = assign_pos_to_range(mids, bin_spans, ref_id="bin_id", pos_col="POS0")
        return mids["bin_id"].fillna(-1).to_numpy(np.int64)

    bin_spans = (
        scaffold.groupby("bin_id", sort=True)
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
    for blk, win_a, dp_a in zip(assay_blocks, window_df_list, dp_corrected_list):
        bin_ids = _windows_to_bins(win_a, bin_spans)
        valid = bin_ids >= 0
        n_drop = int((~valid).sum())
        if n_drop:
            logging.info(
                f"{blk['assay']}: {n_drop}/{len(win_a)} windows outside all bins (dropped)"
            )
        vb = bin_ids[valid]
        win_lengths = (win_a["END"] - win_a["START"]).to_numpy(dtype=np.float64)[valid]
        total_len_per_bin = np.bincount(vb, weights=win_lengths, minlength=num_bbs)
        for s in range(blk["n"]):
            weighted_sums = np.bincount(
                vb, weights=dp_a[valid, s] * win_lengths, minlength=num_bbs
            )
            bb_bases[:, blk["offset"] + s] = weighted_sums
            with np.errstate(invalid="ignore"):
                bb_dp[:, blk["offset"] + s] = weighted_sums / total_len_per_bin
    return bb_dp, bb_bases


def compute_bb_rdr(
    assay_blocks,
    window_df_list,
    dp_corrected_list,
    bb_dp,
    tumor_cols_all,
    base_map,
    rdr_outlier_quantile,
    col_repid,
):
    """Per-bin RDR for every tumor column.

    Each tumor with an RDR base column (from ``base_map``) is normalized by that
    base column, library-size corrected; a tumor without a base is
    median-centered. The base may live in any column (e.g. a different
    assay/platform), so library sizes are computed globally per column. Entries
    above the ``1 - rdr_outlier_quantile`` quantile are set to NaN. Returns a
    ``(num_bbs, len(tumor_cols_all))`` array column-aligned to ``tumor_cols_all``.
    """
    num_bbs, total_samples = bb_dp.shape
    bb_rdr = np.full((num_bbs, len(tumor_cols_all)), np.nan, dtype=np.float32)
    rdr_pos = {c: i for i, c in enumerate(tumor_cols_all)}

    # global per-column total aligned bases for library-size correction
    col_total_bases = np.full(total_samples, np.nan, dtype=np.float64)
    for blk, win_a, dp_a in zip(assay_blocks, window_df_list, dp_corrected_list):
        win_sizes = (win_a["END"] - win_a["START"]).to_numpy(dtype=np.float64)
        tb = np.nansum(dp_a * win_sizes[:, None], axis=0)
        for s in range(blk["n"]):
            col_total_bases[blk["offset"] + s] = tb[s]

    for c in tumor_cols_all:
        m = base_map.get(c)
        if m is not None:
            lib = col_total_bases[m] / col_total_bases[c]
            logging.info(
                f"  bb RDR {col_repid[c]} / base {col_repid[m]}: library factor={lib:.4f}"
            )
            with np.errstate(invalid="ignore", divide="ignore"):
                bb_rdr[:, rdr_pos[c]] = bb_dp[:, c] / bb_dp[:, m] * lib
        else:
            col = bb_dp[:, c]
            valid_i = np.isfinite(col) & (col > 0)
            if valid_i.any():
                med = np.median(col[valid_i])
                logging.info(f"  bb median-centering {col_repid[c]}: median={med:.4f}")
                with np.errstate(invalid="ignore", divide="ignore"):
                    bb_rdr[valid_i, rdr_pos[c]] = col[valid_i] / med

    if rdr_outlier_quantile > 0:
        rdr_upper = np.nanquantile(bb_rdr, 1 - rdr_outlier_quantile)
        n_outlier = int(np.nansum(bb_rdr > rdr_upper))
        logging.info(
            f"RDR outlier filter: quantile={rdr_outlier_quantile}, "
            f"threshold={rdr_upper:.4f}, {n_outlier} entries set to NaN"
        )
        bb_rdr[bb_rdr > rdr_upper] = np.nan
    return bb_rdr


##################################################
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
def assign_snp_bounderies(
    snps: pd.DataFrame, regions: pd.DataFrame, colname="region_id"
):
    """Split each region into one ``[START, END)`` sub-interval per SNP it contains.

    A SNP owns the span between the midpoints of its neighbours, bounded by the
    region edges; a lone SNP owns the whole region. SNPs outside every region keep
    ``START = END = 0`` and an empty ``region_id``. Membership is by 1-based ``POS``
    against the region's coordinates.

    Args:
        snps: SNPs with ``#CHR``, ``POS``, ``POS0``, sorted genomically.
        regions: Intervals with ``#CHR``, ``START``, ``END``, ``region_id`` and
            optionally ``seg_id``.
        colname: Column to write the region identifier into.

    Returns:
        *snps* with ``START``, ``END``, ``BLOCKSIZE``, *colname* (and ``seg_id``).
    """
    has_seg = "seg_id" in regions.columns
    regions = regions.reset_index(drop=True)
    regions["_reg_row"] = np.arange(len(regions))

    snps = assign_pos_to_range(snps, regions, ref_id="_reg_row", pos_col="POS")
    snps["START"] = 0
    snps["END"] = 0
    snps[colname] = ""
    if has_seg:
        snps["seg_id"] = ""

    inside = snps["_reg_row"].notna()
    n_out = int((~inside).sum())
    if n_out:
        logging.info(f"SNP boundaries: {n_out}/{len(snps)} SNPs outside every region")

    reg_start = regions["START"].to_numpy()
    reg_end = regions["END"].to_numpy()
    for reg_row, grp in snps.loc[inside].groupby("_reg_row", sort=False):
        idx = grp.index.to_numpy()
        r = int(reg_row)
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

    snps.drop(columns="_reg_row", inplace=True)
    snps["BLOCKSIZE"] = snps["END"] - snps["START"]
    return snps
