"""GTF feature annotation and per-assay counting into bins.

Last update: 2026-08-11

Functions:
- annotate_feature_type: stamp SNPs with overlapping genes and exon/intron/intergenic
- merge_feature_ids, explode_feature_ids: collapse or expand the joined feature_id
- assign_features_to_ranges: each gene to the range it overlaps most
- read_gene_counts: RNA h5ad UMIs as a (gene, cell) matrix, un-binned
- sum_umis_to_bins: RNA h5ad UMIs into a (bb, cell) matrix
- sum_atac_fragments_to_bins: deduped 10x fragments into a (bb, cell) matrix
"""

import logging

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from io_utils import read_chunks_from_atac_fragments, read_GTF
from segmentation_utils import sum_features_to_bbs
from range_utils import (
    assign_pos_to_range_ovlp,
    assign_range_to_range,
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
    gtf = read_GTF(gtf_file, ("gene", "exon"))
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


def read_gene_counts(h5ad_file, barcodes):
    """Read an RNA h5ad as a gene x cell count matrix, the un-binned RNA unit.

    A gene is indivisible, so it is the RNA counterpart of the window: the finest grid
    ``sum_umis_to_bins`` can aggregate from. Observations are reordered to ``barcodes``
    so the columns match that assay's ``snp.*allele.npz`` matrices.

    Parameters
    ----------
    h5ad_file : str
        AnnData (cells x genes) with ``var`` carrying ``#CHR``, ``START``, ``END``.
    barcodes : sequence of str
        Cell barcodes (``"{raw}_{dataset_id}_{assay_type}"``) in matrix-observation order.

    Returns
    -------
    genes : pd.DataFrame
        ``#CHR``, ``START``, ``END``, ``feature_id`` and ``region_id`` when present, in
        the h5ad's var order.
    x_count : scipy.sparse.csr_matrix
        Shape ``(n_genes, n_cells)``, dtype int32, row-aligned to *genes*.
    """
    import anndata

    adata = anndata.read_h5ad(h5ad_file)
    barcodes = np.asarray(barcodes, dtype=str)
    missing = barcodes[~np.isin(barcodes, adata.obs_names)]
    assert len(missing) == 0, (
        f"h5ad, {len(missing)} barcode(s) missing, e.g. {missing[:5]}"
    )
    adata = adata[barcodes, :].copy()
    genes = adata.var.reset_index(names="feature_id")
    cols = ["#CHR", "START", "END", "feature_id"]
    cols += [c for c in ("region_id",) if c in genes.columns]
    x_count = csr_matrix(adata.X.T).astype(np.int32)
    return genes[cols], x_count


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
        Cell barcodes (``"{raw}_{dataset_id}_{assay_type}"``) in matrix-observation order (that assay's
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
    frag_files, dataset_ids, barcodes, bb_ranges, num_bbs, chunksize=5_000_000
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
    barcodes : pd.DataFrame
        Columns ``dataset_id`` and ``raw`` (the barcode as the fragment file spells it),
        one row per observation and in observation order (identical to that assay's
        ``bb.*allele.npz`` columns), as ``read_barcodes_by_dataset`` returns it.
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
    n_cells = len(barcodes)
    bc_dataset = barcodes["dataset_id"].to_numpy().astype(str)
    bc_raw = barcodes["raw"].to_numpy().astype(str)
    # global observation index keyed by (dataset_id, raw barcode), as the fragments spell it
    obs_of = {(d, r): i for i, (d, r) in enumerate(zip(bc_dataset, bc_raw))}

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
