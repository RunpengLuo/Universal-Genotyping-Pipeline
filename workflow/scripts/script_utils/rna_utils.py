"""RNA / AnnData helpers: gene-to-range assignment and per-bb UMI counts.

The only modules that touch scanpy. A gene is assigned to the range it overlaps most,
never split, so a bb's expression is the sum over whole genes.
"""

import logging

import numpy as np
import pandas as pd
import scanpy as sc

from matrix_utils import sum_features_to_bbs
from range_utils import assign_range_to_range


def assign_features_to_ranges(
    adata: sc.AnnData,
    ranges: pd.DataFrame,
    assay_type: str,
    feature_df_idx="feature_df_idx",
    range_id="region_id",
    drop_cols=True,
):
    """Assign each AnnData feature to the range it overlaps most; drop the rest.

    Features outside every range (masked regions, e.g. centromeres) are removed from
    *adata*, so the returned object carries only assignable genes.

    Args:
        adata: AnnData whose ``var`` carries ``#CHR``, ``START``, ``END``.
        ranges: Target ranges with ``#CHR``, ``START``, ``END`` and *range_id*.
        assay_type: Assay name, for the log lines only.
        feature_df_idx: Temporary column holding each feature's ``var`` row index.
        range_id: Identifier column of *ranges*, carried onto ``var``.
        drop_cols: Drop the two helper columns from ``var`` before returning.

    Returns:
        A copy of *adata* restricted to the assignable features.
    """
    logging.info(f"assign {assay_type} features to ranges, {feature_df_idx}-{range_id}")
    if adata.is_view:
        adata = adata.copy()
    adata.var[feature_df_idx] = np.arange(len(adata.var))

    feature_df = adata.var.reset_index(drop=True)
    logging.info(f"#{assay_type}-features (raw)={len(feature_df)}")

    feature_df = assign_range_to_range(feature_df, ranges, range_id)
    isna_features = feature_df[range_id].isna()
    logging.info(
        f"#{assay_type} feature outside any range={np.sum(isna_features) / len(feature_df):.3%}"
    )
    feature_df.dropna(subset=range_id, inplace=True)
    feature_df[range_id] = feature_df[range_id].astype(ranges[range_id].dtype)
    logging.info(f"#{assay_type} feature (remain)={len(feature_df)}")

    ##################################################
    adata.var = (
        adata.var.reset_index(drop=False)
        .merge(
            right=feature_df[[feature_df_idx, range_id]],
            on=feature_df_idx,
            how="left",
        )
        .set_index("index")
    )
    adata = adata[:, adata.var[range_id].notna()].copy()
    if drop_cols:
        adata.var.drop(columns=[feature_df_idx, range_id], inplace=True)
    else:
        adata.var[range_id] = adata.var[range_id].astype(feature_df[range_id].dtype)
    return adata


def sum_rna_counts_to_bbs(h5ad_file, barcodes, bb_df, num_bbs, assay_type):
    """Aggregate per-cell RNA counts (h5ad from ``process_rna_anndata``) into bbs.

    Each gene is assigned to the bb it overlaps most (``assign_features_to_ranges`` ->
    largest overlap, the same mapping copytyping's ``combine_counts_fixed_bins`` uses);
    its per-cell counts are summed into that bb. Observations are reordered to
    ``barcodes`` so they match that assay's ``bb.*allele.npz`` matrices.

    Parameters
    ----------
    h5ad_file : str
        AnnData (cells x genes) with ``var`` carrying ``#CHR``, ``START``, ``END``.
    barcodes : sequence of str
        Cell barcodes (``"{raw}_{rep}"``) in matrix-observation order (that assay's
        allele observations).
    bb_df : pd.DataFrame
        bbs with ``#CHR``, ``START``, ``END`` (0-based half-open) and ``bb_id``.
    num_bbs : int
        Number of bbs (output features).

    Returns
    -------
    scipy.sparse.csr_matrix, shape ``(num_bbs, n_cells)``, dtype int32.
    """
    adata = sc.read_h5ad(h5ad_file)
    barcodes = np.asarray(barcodes, dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    if len(missing):
        raise ValueError(
            f"{len(missing)} barcodes missing from {h5ad_file}, e.g. {missing[:5]}"
        )
    adata = adata[barcodes, :].copy()
    adata = assign_features_to_ranges(
        adata, bb_df, assay_type, range_id="bb_id", drop_cols=False
    )
    x_count = sum_features_to_bbs(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    return x_count.astype(np.int32)
