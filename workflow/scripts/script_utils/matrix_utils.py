"""Clustered sums over a count matrix, on either axis.

Every matrix in the pipeline is ``(n_features, n_observations)``: a feature is the
genomic entry (bin, bb, SNP, gene), an observation is the measurement entry (cell,
spot, dataset, bulk sample). Binning sums FEATURES (SNPs/bins -> bbs), pseudobulking
sums OBSERVATIONS (cells -> datasets); both are one sparse one-hot multiply, so they
share ``cluster_sum``. A leaf module: numpy and scipy only, importable from anywhere.
"""

import numpy as np
from scipy.sparse import csr_matrix, issparse


def dense_observation(mat, i: int):
    """Observation *i* of a dense or sparse matrix, as a 1-D array over features."""
    obs = mat[:, i]
    return obs.toarray().ravel() if issparse(obs) else np.asarray(obs).ravel()


def cluster_sum(X, cluster_ids, n_clusters, axis=0):
    """Sum the features (``axis=0``) or observations (``axis=1``) of *X* within each cluster.

    Args:
        X: ``(n_features, n_observations)`` dense or sparse matrix.
        cluster_ids: Cluster id per feature (``axis=0``) or per observation
            (``axis=1``), in ``[0, n_clusters)``.
        n_clusters: Number of clusters, i.e. the size of the collapsed axis in the output.
        axis: Axis to collapse.

    Returns:
        ``(n_clusters, n_observations)`` for ``axis=0``, ``(n_features, n_clusters)``
        for ``axis=1``; sparse when *X* is sparse.

    Raises:
        ValueError: *cluster_ids* has the wrong length or holds an id outside the range.
    """
    X = X.tocsr() if issparse(X) else np.asarray(X)
    cluster_ids = np.asarray(cluster_ids, dtype=np.int64)
    n = X.shape[axis]
    if cluster_ids.shape[0] != n:
        raise ValueError(
            f"cluster_ids length {cluster_ids.shape[0]} != axis-{axis} size {n}"
        )
    if n and (cluster_ids.min() < 0 or cluster_ids.max() >= n_clusters):
        raise ValueError("cluster_ids out of range")

    onehot = csr_matrix(
        (np.ones(n, dtype=np.int8), (cluster_ids, np.arange(n, dtype=np.int64))),
        shape=(n_clusters, n),
    )
    return onehot @ X if axis == 0 else X @ onehot.T


def sum_features_to_bbs(X, bb_ids, n_bbs):
    """Sum an SNP- or bin-level matrix into ``(n_bbs, n_observations)``."""
    return cluster_sum(X, bb_ids, n_bbs, axis=0)


def sum_observations_to_pseudobulk(mat, cluster_ids, n_clusters):
    """Sum observations into ``(n_features, n_clusters)`` pseudobulks, always dense."""
    out = cluster_sum(mat, cluster_ids, n_clusters, axis=1)
    return out.toarray() if issparse(out) else np.asarray(out)
