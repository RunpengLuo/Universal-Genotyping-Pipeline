"""RNA / AnnData helpers: gene-to-bin assignment and per-bin UMI counts.

The only modules that touch scanpy. A gene is assigned to the bin it overlaps most,
never split, so a bin's expression is the sum over whole genes.
"""

import logging

import numpy as np
import pandas as pd
import scanpy as sc

from matrix_utils import matrix_segmentation
from interval_utils import assign_interval_to_range


def feature_to_blocks(
    adata: sc.AnnData,
    blocks: pd.DataFrame,
    assay_type: str,
    feature_idx="feature_idx",
    block_idx="region_id",
    drop_cols=True,
):
    """
    filter features not in blocks, likely masked regions include centromeres
    """
    logging.info(f"assign {assay_type} features to blocks, {feature_idx}-{block_idx}")
    if adata.is_view:
        adata = adata.copy()
    adata.var[feature_idx] = np.arange(len(adata.var))

    feature_df = adata.var.reset_index(drop=True)
    logging.info(f"#{assay_type}-features (raw)={len(feature_df)}")

    feature_df = assign_interval_to_range(feature_df, blocks, block_idx)
    isna_features = feature_df[block_idx].isna()
    logging.info(
        f"#{assay_type} feature outside any blocks={np.sum(isna_features) / len(feature_df):.3%}"
    )
    feature_df.dropna(subset=block_idx, inplace=True)
    feature_df[block_idx] = feature_df[block_idx].astype(blocks[block_idx].dtype)
    logging.info(f"#{assay_type} feature (remain)={len(feature_df)}")

    ##################################################
    adata.var = (
        adata.var.reset_index(drop=False)
        .merge(
            right=feature_df[[feature_idx, block_idx]],
            on=feature_idx,
            how="left",
        )
        .set_index("index")
    )
    adata = adata[:, adata.var[block_idx].notna()].copy()
    if drop_cols:
        adata.var.drop(columns=[feature_idx, block_idx], inplace=True)
    else:
        adata.var[block_idx] = adata.var[block_idx].astype(feature_df[block_idx].dtype)
    return adata


def rna_h5ad_to_bb(h5ad_file, barcodes, bb_df, num_bbs, assay_type):
    """Aggregate per-cell RNA counts (h5ad from ``process_rna_anndata``) into bb bins.

    Each RNA feature (gene) is assigned to the bb bin it overlaps most (``feature_to_blocks``
    -> largest overlap, same mapping as copytyping's ``combine_counts_fixed_bins``); its per-cell counts
    are summed into that bin. Cells are reordered to ``barcodes`` so the columns match that
    assay's ``bb.*allele.npz`` matrices.

    Parameters
    ----------
    h5ad_file : str
        AnnData (cells x genes) with ``var`` carrying ``#CHR``, ``START``, ``END``.
    barcodes : sequence of str
        Cell barcodes (``"{raw}_{rep}"``) in matrix-column order (that assay's allele columns).
    bb_df : pd.DataFrame
        Bins with ``#CHR``, ``START``, ``END`` (0-based half-open) and ``bb_id``.
    num_bbs : int
        Number of bins (output rows).

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
    adata = feature_to_blocks(
        adata, bb_df, assay_type, block_idx="bb_id", drop_cols=False
    )
    x_count = matrix_segmentation(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    return x_count.astype(np.int32)
