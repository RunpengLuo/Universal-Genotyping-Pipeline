"""GTF feature annotation and per-assay feature counting.

A feature here is a GTF entity, a gene or an exon. Two groups:

- annotation - stamp SNPs with the genes they sit in (``annotate_feature_type``), collapse
  ``feature_id`` strings (``merge_feature_ids``), and glue each gene's span of fixed bins
  into one cluster so no bb splits a gene (``stamp_gene_clusters``). The filter on the
  result is ``phase_and_concat_utils.get_mask_by_exon``, next to the other SNP masks.
- counting - turn an assay's raw records into a ``(bin, cell)`` count matrix:
  ``sum_umis_to_bins`` for the scRNA/VISIUM h5ad, ``sum_atac_fragments_to_bins`` for 10x
  fragment files. ``assign_features_to_ranges`` is the gene-to-range mapping both the RNA
  path and ``combine_counts_fixed_bins`` use; a gene is assigned to the range it overlaps
  most and never split, so a bb's expression is the sum over whole genes.

anndata is imported inside the one function that reads an h5ad, so the module imports
without the single-cell stack.
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from io_utils import read_chunks_from_atac_fragments, read_gtf
from matrix_utils import sum_features_to_bbs
from range_utils import (
    assign_pos_to_range_ovlp,
    assign_range_to_range,
    merge_ranges_to_clusters,
    overlaps_any_range,
)
from utils import add_chr_prefix


def merge_feature_ids(strings, sep=";", default="intergenic"):
    """Collapse an iterable of ``sep``-joined feature_id strings into one deduped union.

    Drops ``default`` tokens unless nothing else remains; preserves first-seen order.
    """
    seen = dict()
    for s in strings:
        if not isinstance(s, str):
            continue
        for tok in s.split(sep):
            if tok and tok != default:
                seen[tok] = None
    return sep.join(seen) if seen else default


def annotate_feature_type(snps, gtf_file):
    """Annotate SNPs with feature_id (overlapping genes) and feature_type.

    ``feature_id`` is a ``;``-joined list of every GTF gene the SNP overlaps
    (``intergenic`` when none), so gene membership is read back off it rather than
    assigned a second time. ``feature_type`` is exon > intron > intergenic.
    """
    gtf = read_gtf(gtf_file, ("gene", "exon"))
    snps, _ = assign_pos_to_range_ovlp(
        snps, gtf["gene"], ref_id="gene_id", out_col="feature_id", fillna="intergenic"
    )
    in_gene = snps["feature_id"] != "intergenic"
    in_exon = pd.Series(overlaps_any_range(snps, gtf["exon"]), index=snps.index)

    snps["feature_type"] = "intergenic"
    snps.loc[in_gene, "feature_type"] = "intron"
    snps.loc[in_exon, "feature_type"] = "exon"
    return snps


def explode_feature_ids(df, cols=None, sep=";"):
    """One row per (row, gene) from the *sep*-joined ``feature_id`` column.

    Rows with no gene (NaN or ``intergenic``) are dropped, so the result carries only
    real gene ids.

    Args:
        df: Frame with a ``feature_id`` column.
        cols: Columns to keep besides ``feature_id``; ``None`` keeps all.
        sep: Separator joining the gene ids, matching ``merge_feature_ids``.

    Returns:
        A copy of the kept rows, one per gene.
    """
    genic = df[df["feature_id"].notna() & (df["feature_id"] != "intergenic")]
    genic = genic[list(cols) + ["feature_id"]].copy() if cols else genic.copy()
    genic["feature_id"] = genic["feature_id"].str.split(sep)
    genic = genic.explode("feature_id")
    return genic[genic["feature_id"] != "intergenic"]


def stamp_gene_clusters(bin_df, snps_binned):
    """Glue each gene's span of fixed bins into one cluster, so no bb splits a gene.

    The ``;``-joined multi-gene ``feature_id`` is exploded so each gene gets its own
    first..last bin span. Spans cover the SNP-free bins between a gene's SNPs too, so a
    bb boundary cannot land in an intronic gap.

    Args:
        bin_df: Fixed bins with ``bin_id``. Modified in place.
        snps_binned: SNPs carrying ``bin_id`` and ``feature_id``.

    Returns:
        *bin_df* with ``gene_cluster`` added.
    """
    genic = explode_feature_ids(snps_binned, cols=["bin_id"])
    rng = genic.groupby("feature_id")["bin_id"].agg(["min", "max"])
    bin_df["gene_cluster"] = merge_ranges_to_clusters(
        len(bin_df), zip(rng["min"].to_numpy(), rng["max"].to_numpy() + 1)
    )
    logging.info(
        f"gene-aware binning: {len(rng)} genes over {len(bin_df)} fixed bins -> "
        f"{bin_df['gene_cluster'].nunique()} gene/intergenic clusters (bbs never split a gene)"
    )
    return bin_df


def assign_features_to_ranges(
    adata,
    ranges: pd.DataFrame,
    assay_type: str,
    range_id="region_id",
):
    """Assign each AnnData feature to the range it overlaps most; drop the rest.

    Features outside every range (masked regions, e.g. centromeres) are removed from
    *adata*, so the returned object carries only assignable genes.

    Args:
        adata: AnnData whose ``var`` carries ``#CHR``, ``START``, ``END``.
        ranges: Target ranges with ``#CHR``, ``START``, ``END`` and *range_id*.
        assay_type: Assay name, for the log lines only.
        range_id: Identifier column of *ranges*, carried onto ``var``.

    Returns:
        A copy of *adata* restricted to the assignable features, ``var`` carrying
        *range_id* in *ranges*' dtype.
    """
    logging.info(f"assign {assay_type} features to ranges, {range_id}")
    if adata.is_view:
        adata = adata.copy()
    n_raw = len(adata.var)
    logging.info(f"#{assay_type}-features (raw)={n_raw}")

    # var rows and the returned rows are in the same order, so assign positionally
    annotated, na_idx = assign_range_to_range(adata.var, ranges, range_id)
    adata.var[range_id] = annotated[range_id].to_numpy()
    logging.info(
        f"#{assay_type} feature outside any range={len(na_idx) / max(n_raw, 1):.3%}"
    )

    adata = adata[:, adata.var[range_id].notna()].copy()
    adata.var[range_id] = adata.var[range_id].astype(ranges[range_id].dtype)
    logging.info(f"#{assay_type} feature (remain)={adata.n_vars}")
    return adata


def sum_umis_to_bins(h5ad_file, barcodes, bb_df, num_bbs, assay_type):
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
        Cell barcodes (``"{raw}_{dataset_id}"``) in matrix-observation order (that assay's
        allele observations).
    bb_df : pd.DataFrame
        bbs with ``#CHR``, ``START``, ``END`` (0-based half-open) and ``bb_id``.
    num_bbs : int
        Number of bbs (output features).

    Returns
    -------
    scipy.sparse.csr_matrix, shape ``(num_bbs, n_cells)``, dtype int32.
    """
    import anndata

    adata = anndata.read_h5ad(h5ad_file)
    barcodes = np.asarray(barcodes, dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    assert len(missing) == 0, (
        f"h5ad, {len(missing)} barcode(s) missing, e.g. {missing[:5]}"
    )
    adata = adata[barcodes, :].copy()
    adata = assign_features_to_ranges(adata, bb_df, assay_type, range_id="bb_id")
    x_count = sum_features_to_bbs(adata.X.T, adata.var["bb_id"].to_numpy(), num_bbs)
    return x_count.astype(np.int32)


def sum_atac_fragments_to_bins(
    frag_files, dataset_ids, barcodes_full, bb_ranges, num_bbs, chunksize=5_000_000
):
    """Count deduped ATAC fragments per bb per cell from 10x fragment files.

    Each record of a 10x ``atac_fragments.tsv.gz`` is one deduplicated fragment
    (``chrom, start, end, barcode, readSupport``); the readSupport field is IGNORED.
    Every fragment is counted once, assigned to the bb containing its midpoint, so
    each observation sums to that cell's in-bb fragment count.

    Parameters
    ----------
    frag_files, dataset_ids : parallel lists
        ``frag_files[i]`` is the fragment file for replicate ``dataset_ids[i]``.
    barcodes_full : pd.DataFrame
        Columns ``REP_ID``, ``BARCODE`` (``BARCODE`` = ``"{raw}_{dataset_id}"``) giving the observation
        order (identical to that assay's ``bb.*allele.npz`` observations).
    bb_ranges : pd.DataFrame
        Non-overlapping ``#CHR``, ``START``, ``END`` (0-based half-open) ranges carrying
        ``bb_id``; several ranges may share a ``bb_id``. Pass the FIXED BINS stamped with
        their owning bb rather than the bb hulls: a fragment in a blacklist hole then hits
        no bin and is dropped, where a hull would swallow it.
    num_bbs : int
        Number of bbs (output features).

    Returns
    -------
    scipy.sparse.csr_matrix, shape ``(num_bbs, n_cells)``, dtype int32.
    """
    n_cells = len(barcodes_full)
    bc_dataset = barcodes_full["REP_ID"].to_numpy().astype(str)
    bc_full = barcodes_full["BARCODE"].to_numpy().astype(str)
    # global observation index keyed by (dataset_id, raw_barcode); strip the "_{dataset_id}" suffix
    obs_of = {}
    for i in range(n_cells):
        dataset_id, raw = bc_dataset[i], bc_full[i]
        sfx = "_" + dataset_id
        if raw.endswith(sfx):
            raw = raw[: -len(sfx)]
        obs_of[(dataset_id, raw)] = i

    bb_all, obs_all = [], []
    for frag_file, dataset_id in zip(frag_files, dataset_ids):
        dataset_map = {raw: c for (r, raw), c in obs_of.items() if r == dataset_id}
        if not dataset_map or frag_file is None:
            continue
        n_frag = n_seen = 0
        for chunk in read_chunks_from_atac_fragments(frag_file, chunksize=chunksize):
            obs_vals = chunk["BC"].map(dataset_map).to_numpy()
            m = ~pd.isna(obs_vals)
            if not m.any():
                continue
            sub = chunk.loc[m]
            frag = pd.DataFrame(
                {
                    "#CHR": add_chr_prefix(sub["#CHR"]).to_numpy(),
                    "START": sub["start"].to_numpy(),
                    "END": sub["end"].to_numpy(),
                }
            )
            frag, _ = assign_range_to_range(frag, bb_ranges, "bb_id", rule="midpoint")
            keep = frag["bb_id"].notna().to_numpy()
            if not keep.any():
                continue
            bb_all.append(frag.loc[keep, "bb_id"].to_numpy().astype(np.int64))
            obs_all.append(obs_vals[m][keep].astype(np.int64))
            n_frag += int(keep.sum())
            n_seen += int(m.sum())
        pct = 100.0 * n_frag / n_seen if n_seen else 0.0
        logging.info(
            f"  ATAC {dataset_id}: {n_frag}/{n_seen} ({pct:.1f}%) fragments counted; "
            "the rest fall outside every bin (blacklist holes, off-segment)"
        )

    if bb_all:
        bbs = np.concatenate(bb_all)
        obs = np.concatenate(obs_all)
    else:
        bbs = np.zeros(0, dtype=np.int64)
        obs = np.zeros(0, dtype=np.int64)
    data = np.ones(len(bbs), dtype=np.int32)
    return csr_matrix((data, (bbs, obs)), shape=(num_bbs, n_cells), dtype=np.int32)
