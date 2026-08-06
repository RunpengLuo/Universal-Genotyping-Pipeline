"""Grouped sums over a count matrix, on either axis.

Binning sums ROWS (SNPs/windows -> bins), pseudobulking sums COLUMNS
(cells -> replicates); both are one sparse one-hot multiply, so they share
``group_sum``. A leaf module: numpy and scipy only, importable from anywhere.
"""

import numpy as np
from scipy.sparse import csr_matrix, issparse


def dense_col(mat, i: int):
    """Column *i* of a dense or sparse matrix, as a 1-D array."""
    col = mat[:, i]
    return col.toarray().ravel() if issparse(col) else np.asarray(col).ravel()


def group_sum(X, ids, n_groups, axis=0):
    """Sum the rows (``axis=0``) or columns (``axis=1``) of *X* within each group.

    Args:
        X: ``(N, M)`` dense or sparse matrix.
        ids: Group id per row (``axis=0``) or per column (``axis=1``), in ``[0, n_groups)``.
        n_groups: Number of groups, i.e. the size of the grouped axis in the output.
        axis: Axis to collapse.

    Returns:
        ``(n_groups, M)`` for ``axis=0``, ``(N, n_groups)`` for ``axis=1``; sparse
        when *X* is sparse.

    Raises:
        ValueError: *ids* has the wrong length or holds an id outside the range.
    """
    X = X.tocsr() if issparse(X) else np.asarray(X)
    ids = np.asarray(ids, dtype=np.int64)
    n = X.shape[axis]
    if ids.shape[0] != n:
        raise ValueError(f"ids length {ids.shape[0]} != axis-{axis} size {n}")
    if n and (ids.min() < 0 or ids.max() >= n_groups):
        raise ValueError("ids out of range")

    onehot = csr_matrix(
        (np.ones(n, dtype=np.int8), (ids, np.arange(n, dtype=np.int64))),
        shape=(n_groups, n),
    )
    return onehot @ X if axis == 0 else X @ onehot.T


def matrix_segmentation(X, bin_ids, K):
    """Sum an ``(N, M)`` SNP/window-by-sample matrix into ``(K, M)`` bins."""
    return group_sum(X, bin_ids, K, axis=0)


def pseudobulk_by_groups(mat, group_idx, n_groups):
    """Sum a feature-by-cell matrix into ``(n_features, n_groups)``, always dense."""
    out = group_sum(mat, group_idx, n_groups, axis=1)
    return out.toarray() if issparse(out) else np.asarray(out)
