"""Read the pipeline's input files and write its tabular outputs.

Last update: 2026-08-12

Functions:
- read_VCF, read_BED, read_segment_bed: parse the coordinate inputs
- read_gtf, read_genes_gtf_file: parse the gene annotation GTF
- read_chrom_sizes, read_bedgraph, read_window_bed: parse the reference and bin inputs
- read_allele_mat, read_snp_mats: allele matrices, dense bulk or sparse single-cell
- read_barcodes, read_barcodes_by_dataset: the single-cell column axis
- read_chunks_from_atac_fragments: stream a 10x fragment file in chunks
- read_10x_ranger_scRNA, read_10x_ranger_spatial: Cell and Space Ranger matrices
- read_bcftools_pileup_counts: a bcftools AD table as count matrices
- write_snp_info, write_bb_file, write_sample_ids: the three output schemas
"""

import logging
import os
import tempfile
from collections import OrderedDict

import pandas as pd
import numpy as np

from const import GTF_COLUMNS, RANGER_MATRIX_H5, RANGER_SPATIAL_DIR
from utils import add_chr_prefix, sort_chroms, sort_df_chr


def read_chrom_sizes(sz_file: str):
    """Read a two-column chromosome-sizes file.

    Args:
        sz_file: Tab-separated ``chrom<TAB>size``.

    Returns:
        OrderedDict of chromosome name to length, in file order.
    """
    chr_sizes = OrderedDict()
    with open(sz_file, "r") as rfd:
        for line in rfd.readlines():
            ch, sizes = line.strip().split()
            chr_sizes[ch] = int(sizes)
        rfd.close()
    return chr_sizes


def read_VCF(
    vcf_file: str,
    addchr=True,
    addkey=False,
    snps_presorted=False,
    add_pos0=False,
    add_phase1=False,
):
    """Read a VCF into a DataFrame, exploding its INFO and FORMAT fields into columns.

    Args:
        vcf_file: Path to the (optionally gzipped) VCF.
        addchr: Prepend ``chr`` to contigs named without it.
        addkey: Add ``KEY`` (``#CHROM_POS``).
        snps_presorted: Skip the genomic sort.
        add_pos0: Add ``POS0`` (0-based).
        add_phase1: Add ``PHASE``, the second GT allele.

    Returns:
        DataFrame with the 8 fixed VCF columns, ``#CHR``, ``RAW_SNP_DF_IDX`` and one
        column per INFO/FORMAT key; None when the file has no records.
    """
    snps = pd.read_csv(
        vcf_file, comment="#", sep="\t", header=None, dtype={0: "string"}
    )
    if snps.empty:
        return None
    ncols = snps.shape[1]
    assert ncols == 8 or ncols >= 10, (
        f"VCF file, expected 8 or >=10 columns, got {ncols}"
    )
    colnames = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if ncols >= 10:
        colnames += ["FORMAT", "SAMPLE"]
        snps = snps.iloc[:, :10].copy()
    snps.columns = colnames
    snps["POS"] = snps["POS"].astype(np.int64)
    snps["RAW_SNP_DF_IDX"] = np.arange(len(snps))
    if addchr and not str(snps["#CHROM"].iloc[0]).startswith("chr"):
        snps["#CHROM"] = "chr" + snps["#CHROM"].astype(str)

    snps["#CHROM"] = snps["#CHROM"].str.replace("^chrMT$", "chrM", regex=True)
    chrom_order = sort_chroms(snps["#CHROM"].unique().tolist())
    snps["#CHROM"] = pd.Categorical(
        snps["#CHROM"], categories=chrom_order, ordered=True
    )
    if not snps_presorted:
        snps = snps.sort_values(["#CHROM", "POS"], kind="mergesort")

    snps["#CHR"] = snps["#CHROM"]

    # parse INFO column
    info_kvs = (
        snps["INFO"]
        .fillna("")
        .str.split(";")
        .explode()
        .loc[lambda s: s.ne("")]
        .to_frame("kv")
    )
    kv = info_kvs["kv"].str.split("=", n=1, expand=True)
    info_kvs["key"] = kv[0]
    info_kvs["val"] = kv[1] if kv.shape[1] > 1 else None
    info_kvs["val"] = info_kvs["val"].fillna(True)  # INFO flags (no '=') -> True
    info_kvs["row"] = info_kvs.index

    info_wide = info_kvs.pivot_table(
        index="row", columns="key", values="val", aggfunc="first"
    )
    snps = snps.join(info_wide)

    # parse FORMAT column
    if "FORMAT" in snps.columns:
        fmt_keys = snps["FORMAT"].fillna("").str.split(":")
        samp_vals = snps["SAMPLE"].fillna("").str.split(":")

        fmt_long = pd.DataFrame(
            {"key": fmt_keys.explode(), "val": samp_vals.explode()}
        ).dropna(subset=["key"])
        fmt_long["row"] = fmt_keys.explode().index  # original variant row index

        fmt_wide = fmt_long.pivot_table(
            index="row", columns="key", values="val", aggfunc="first"
        )
        snps = snps.drop(columns=["FORMAT", "SAMPLE"])
        # A key can appear in both INFO and FORMAT (e.g. AD/DP when a phaser
        # such as longphase preserves the input INFO and FORMAT fields). Detect
        # any such overlap and let the per-sample FORMAT value take precedence,
        # dropping the INFO-derived duplicate to avoid a "columns overlap" join
        # error. No overlap (eagle/shapeit emit only GT) -> no-op.
        dup_cols = snps.columns.intersection(fmt_wide.columns)
        if len(dup_cols):
            snps = snps.drop(columns=dup_cols)
        snps = snps.join(fmt_wide)

    if addkey:
        snps["KEY"] = snps["#CHROM"].astype(str) + "_" + snps["POS"].astype(str)
    if add_pos0:
        snps["POS0"] = snps["POS"] - 1
    if add_phase1:
        snps["PHASE"] = snps["GT"].str[2].astype(np.int8).to_numpy()
    snps = snps.reset_index(drop=True)
    return snps


def read_bcftools_pileup_counts(tsv_file: str, parent_alt_by_key: dict):
    """Read a bcftools per-locus AD table as pseudobulk depth/alt count matrices.

    Counts the parent ALT allele only, looked up in the locus's own ALT list, so
    ``DP = ref + alt`` and a downstream ``REF = DP - ALT`` is exact.

    Args:
        tsv_file: Path to the (optionally gzipped) counts TSV.
        parent_alt_by_key: Parent ALT allele keyed by ``#CHROM_POS``.

    Returns:
        ``(snps, tot_mtx, ad_mtx)``: *snps* carries ``KEY`` (matching ``read_VCF``) and
        ``RAW_SNP_DF_IDX`` (file row order); the matrices are ``(len(snps), 1)`` csr,
        shaped for ``map_allele_mat_to_snps``.
    """
    from scipy.sparse import csr_matrix  # scipy is not a runner-env dependency

    df = pd.read_csv(
        tsv_file,
        sep="\t",
        header=None,
        names=["#CHROM", "POS", "REF", "ALT", "AD"],
        dtype={"#CHROM": "string", "REF": "string", "ALT": "string", "AD": "string"},
    )
    chrom = add_chr_prefix(df["#CHROM"]).str.replace("^chrMT$", "chrM", regex=True)
    keys = chrom + "_" + df["POS"].astype(np.int64).astype(str)
    alt_lists = df["ALT"].str.split(",")
    ad_lists = df["AD"].str.split(",").apply(lambda xs: [int(x) for x in xs])
    parent_alts = keys.map(parent_alt_by_key)

    ref = np.array([ad[0] if ad else 0 for ad in ad_lists], dtype=np.int64)
    alt = np.array(
        [
            ad[1 + alts.index(pa)] if pa in alts else 0
            for ad, alts, pa in zip(ad_lists, alt_lists, parent_alts)
        ],
        dtype=np.int64,
    )
    snps = pd.DataFrame({"KEY": keys.to_numpy(), "RAW_SNP_DF_IDX": np.arange(len(df))})
    # a csr built from a dense column drops the zeros itself
    tot_mtx = csr_matrix((ref + alt).reshape(-1, 1))
    ad_mtx = csr_matrix(alt.reshape(-1, 1))
    return snps, tot_mtx, ad_mtx


def read_allele_mat(npz_file, mat_dtype=None):
    """Read one SNP-level allele matrix, dense or sparse, as it was written.

    Bulk writes a dense ``mat`` key, single-cell a ``scipy.sparse`` archive; the keys
    tell them apart.

    Args:
        npz_file: Path to the ``.npz``.
        mat_dtype: Cast the matrix to this dtype; ``None`` keeps the stored one.

    Returns:
        ``np.ndarray`` for a dense file, ``scipy.sparse.csr_matrix`` for a sparse one.
    """
    from scipy.sparse import load_npz  # scipy is not a runner-env dependency

    with np.load(npz_file) as npz:
        is_dense = "mat" in npz.files
        mat = npz["mat"] if is_dense else load_npz(npz_file)
    return mat if mat_dtype is None else mat.astype(mat_dtype)


def read_snp_mats(snp_info_file, tot_file, a_file, b_file, mat_dtype=None):
    """Read a SNP table and its T/A/B allele matrices, in file row order.

    Args:
        snp_info_file: ``snps.tsv.gz`` from phase_and_concat.
        tot_file, a_file, b_file: the total / A-allele / B-allele ``.npz``.
        mat_dtype: Cast the three matrices to this dtype; ``None`` keeps the stored one,
            which for a sparse file avoids copying its ``data`` array.

    Returns:
        ``(snps, tot_mtx, a_mtx, b_mtx)``; the matrices are dense or sparse per
        ``read_allele_mat``.
    """
    snps = pd.read_table(snp_info_file, sep="\t")
    return (
        snps,
        read_allele_mat(tot_file, mat_dtype),
        read_allele_mat(a_file, mat_dtype),
        read_allele_mat(b_file, mat_dtype),
    )


def read_BED(bed_file: str, addchr=True, col_id="region_id"):
    """Read BED file.

    Args:
        bed_file: Path to a BED file, at least 3 columns.
        addchr: Prepend ``chr`` to contigs if contigs are not chr-prefix.
        col_id: Name of the id column, read from column 4 or derived.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END`` and *col_id*.

    Notes/References:
        Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    """
    df = pd.read_table(bed_file, sep="\t", header=None, dtype={0: "string"})
    assert len(df.columns) >= 3, (
        f"{bed_file}: BED file, expected >=3 columns, got {len(df.columns)}"
    )
    has_name = len(df.columns) >= 4
    df = df.iloc[:, : 4 if has_name else 3].copy()
    df.columns = ["#CHR", "START", "END"] + ([col_id] if has_name else [])
    if not str(df["#CHR"].iloc[0]).startswith("chr") and addchr:
        df["#CHR"] = "chr" + df["#CHR"].astype(str)
    if not has_name:
        df[col_id] = (
            df["#CHR"] + ":" + df["START"].astype(str) + "-" + df["END"].astype(str)
        )
    return df


def read_segment_bed(bed_file: str, addchr=True):
    """Read custom BED file with 5th column defines segment ID.

    Args:
        bed_file: Path to the custom BED, 5 columns.
        addchr: Prepend ``chr`` to contigs named without it.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END``, ``region_id``, ``seg_id``.
    """
    df = read_BED(bed_file, addchr=addchr, col_id="region_id")
    df["seg_id"] = pd.read_table(
        bed_file, sep="\t", header=None, usecols=[4], dtype="string"
    ).iloc[:, 0]
    return df


def read_window_bed(bed_files, chroms=None):
    """Read fixed-bin BED(s) into one genomically sorted frame with ``bin_id``.

    Takes one path, or several to union: the bulk path has one ``window.tsv.gz`` per
    assay, each the same tiling minus that assay's NaN bins, so they are unioned on the
    coordinates. The bias-correction covariates (GC/MAP/REPLI) are dropped - binning
    groups by cluster and would copy them per group. ``seg_id`` is carried only when
    every file has it, else it falls back to ``region_id`` (one segment per arm).

    Every producer writes a sorted file, but a union of differing subsets is not sorted,
    hence the sort here. ``bin_id`` must be the row position: ``build_adaptive_bins``
    maps ``bb_id`` back positionally.

    Args:
        bed_files: One path, or a sequence of paths to union.
        chroms: Keep only these contigs; ``None`` keeps every row. Single-cell reads the
            window BED itself, which may be a genome-wide grid, so it filters; bulk
            passes ``None`` because ``rd_correct`` already trimmed each window.tsv.gz.

    Returns:
        ``(bins, raw_bins)``: *bins* is the union, with ``#CHR``, ``START``, ``END``,
        ``region_id``, ``seg_id``, ``bin_id``. *raw_bins* holds each file exactly as
        read, so a caller whose per-file matrix is row-aligned to it (bulk's corrected
        depth) can use it without reading the files a second time.
    """
    if isinstance(bed_files, (str, os.PathLike)):
        bed_files = [bed_files]
    raw_bins = [pd.read_table(f, sep="\t", dtype={"#CHR": str}) for f in bed_files]
    has_seg = all("seg_id" in f.columns for f in raw_bins)
    cols = ["#CHR", "START", "END", "region_id"] + (["seg_id"] if has_seg else [])
    bin_df = pd.concat([f[cols] for f in raw_bins], ignore_index=True)
    bin_df["#CHR"] = add_chr_prefix(bin_df["#CHR"])
    if chroms is not None:
        bin_df = bin_df[bin_df["#CHR"].isin(chroms)]
    bin_df = bin_df.drop_duplicates(["#CHR", "START", "END"])
    bin_df = sort_df_chr(bin_df, ch="#CHR", pos="START").reset_index(drop=True)
    if not has_seg:
        bin_df["seg_id"] = bin_df["region_id"]
    bin_df["bin_id"] = np.arange(len(bin_df))
    return bin_df, raw_bins


def read_bedgraph(bg_file: str, chroms=None):
    """Read a bedGraph track: ``chrom start end value``, 0-based half-open.

    Args:
        bg_file: Path to the (optionally gzipped) bedGraph.
        chroms: Keep only these contigs; ``None`` keeps every row.

    Returns:
        DataFrame with ``#CHR``, ``START``, ``END``, ``signal``, reindexed from 0.

    Notes/References:
        Format: https://genome.ucsc.edu/goldenPath/help/bedgraph.html
    """
    df = pd.read_csv(
        bg_file,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "signal"],
        dtype={"#CHR": str, "START": np.int64, "END": np.int64, "signal": np.float64},
    )
    df["#CHR"] = add_chr_prefix(df["#CHR"])
    if chroms is not None:
        df = df[df["#CHR"].isin(chroms)]
    return df.reset_index(drop=True)


def read_barcodes(bc_file: str):
    """Read a barcode file, one barcode per line.

    Args:
        bc_file: Path to the barcode list.

    Returns:
        List of barcodes, in file order.
    """
    barcodes = (
        pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()
    )
    return barcodes


def read_barcodes_by_dataset(bc_file: str):
    """Read ``barcodes.tsv.gz`` and split each key back into its three fields.

    The parse is positional: no assay_type holds ``_`` and no raw barcode does either
    (asserted at write), so the last and first ``_`` bound the dataset_id.

    Args:
        bc_file: One ``{raw}_{dataset_id}_{assay_type}`` per line, no header, in
            matrix-column order.

    Returns:
        DataFrame with ``raw``, ``dataset_id``, ``assay_type`` and ``BARCODE`` (the full
        key), one row per matrix column, in file order.
    """
    full = pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0]
    head, assay_type = _split_once(full, bc_file, from_right=True)
    raw, dataset_id = _split_once(head, bc_file, from_right=False)
    return pd.DataFrame(
        {
            "raw": raw,
            "dataset_id": dataset_id,
            "assay_type": assay_type,
            "BARCODE": full,
        }
    )


def _split_once(values: pd.Series, bc_file: str, from_right: bool):
    """Split every barcode field on one ``_``, asserting each side is non-empty."""
    parts = (
        values.str.rsplit("_", n=1, expand=True)
        if from_right
        else values.str.split("_", n=1, expand=True)
    )
    assert parts.shape[1] == 2 and parts.notna().all().all(), (
        f"{bc_file}, every barcode must be '{{raw}}_{{dataset_id}}_{{assay_type}}'"
    )
    return parts[0], parts[1]


def read_chunks_from_atac_fragments(frag_file: str, chunksize=5_000_000):
    """Read a 10x ATAC fragment file in chunks, keeping its first four columns.

    Contig names are left as the file spells them: one sample runs to hundreds of
    millions of records, so a caller chr-normalizes the rows it keeps.

    Args:
        frag_file: Path to the (optionally gzipped) fragment TSV.
        chunksize: Records per chunk.

    Returns:
        Iterator of DataFrames with ``#CHR``, ``start``, ``end``, ``BC``.

    Notes/References:
        Format: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/fragments-file
    """
    return pd.read_csv(
        frag_file,
        sep="\t",
        comment="#",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["#CHR", "start", "end", "BC"],
        dtype={0: str, 1: np.int64, 2: np.int64, 3: str},
        chunksize=chunksize,
    )


def read_10x_ranger_spatial(
    matrix_h5, names, paths, library_id, assay_type, load_images=True
):
    """Read one Space Ranger spatial dataset into an AnnData.

    squidpy takes a directory, so the Space Ranger layout is rebuilt as symlinks in a
    temporary directory for the read; names come from ``RANGER_*`` in const.py.

    Args:
        matrix_h5: Path to this dataset's feature-barcode matrix.
        names: Space Ranger filenames under spatial/, for this dataset.
        paths: Paths supplying those files, in the same order.
        library_id: Library id squidpy records in ``uns``; the dataset_id.
        assay_type: VISIUM | VISIUM3prime.
        load_images: Read the tissue images; VISIUM3prime must pass False.

    Returns:
        AnnData with unique var_names.

    Raises:
        AssertionError: load_images is set for VISIUM3prime.

    Notes/References:
        spatial/ layout: https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/outputs/spatial-outputs
    """
    import squidpy as sq

    if assay_type == "VISIUM3prime":
        assert not load_images, "VISIUM3prime, squidpy cannot load its tissue images"

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.symlink(
            os.path.abspath(matrix_h5),
            os.path.join(tmp_dir, RANGER_MATRIX_H5[0]),
        )
        spatial_dir = os.path.join(tmp_dir, RANGER_SPATIAL_DIR)
        os.makedirs(spatial_dir)
        for name, path in zip(names, paths):
            os.symlink(os.path.abspath(path), os.path.join(spatial_dir, name))
        logging.info(f"staged {len(names) + 1} files for squidpy: {names}")
        adata = sq.read.visium(tmp_dir, load_images=load_images, library_id=library_id)
    adata.var_names_make_unique()
    return adata


def read_10x_ranger_scRNA(matrix_h5):
    """Read one Cell Ranger gene-expression matrix into an AnnData.

    Args:
        matrix_h5: Path to ``filtered_feature_bc_matrix.h5``.

    Returns:
        AnnData of the gene-expression features only, with unique var_names.

    Notes/References:
        Format: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-h5-matrices
    """
    import scanpy as sc

    adata = sc.read_10x_h5(matrix_h5, gex_only=True)
    adata.var_names_make_unique()
    return adata


def read_gtf(gtf_file: str, feature_types):
    """Parse a GTF once and split it by feature type.

    Contigs are chr-normalized and coordinates become 0-based half-open; genes are
    deduplicated by ``gene_id``.

    Args:
        gtf_file: Path to a GTF annotation file (optionally gzipped).
        feature_types: Feature types to extract, e.g. ``("gene", "exon")``.

    Returns:
        ``{feature_type: DataFrame}`` with ``#CHR``, ``START``, ``END``, ``gene_id``.

    Notes/References:
        Format: https://genome.ucsc.edu/FAQ/FAQformat.html#format4
    """
    wanted = list(feature_types)
    gtf = pd.read_csv(
        gtf_file,
        sep="\t",
        comment="#",
        header=None,
        names=GTF_COLUMNS,
        dtype={"seqname": str},
        low_memory=False,
    )
    gtf = gtf.loc[
        gtf["feature"].isin(wanted),
        ["feature", "seqname", "start", "end", "attributes"],
    ]
    flat = pd.DataFrame(
        {
            "feature": gtf["feature"].values,
            "#CHR": add_chr_prefix(gtf["seqname"]).values,
            "START": gtf["start"].values - 1,  # GTF is 1-based -> 0-based
            "END": gtf["end"].values,  # GTF end is inclusive -> half-open
            "gene_id": gtf["attributes"]
            .str.extract(r'gene_id "([^"]+)"', expand=False)
            .values,
        }
    )
    out = {}
    for feature_type in wanted:
        sub = flat.loc[flat["feature"] == feature_type].drop(columns="feature")
        if feature_type == "gene":
            sub = sub.drop_duplicates("gene_id", keep="first")
        out[feature_type] = sub.reset_index(drop=True)
    return out


def read_genes_gtf_file(gtf_file: str, id_col="gene_ids"):
    """Gene-level GTF records with 0-based coordinates, deduplicated by gene.

    Args:
        gtf_file: Path to a GTF annotation file.
        id_col: Column name for the gene identifier in the output.

    Returns:
        DataFrame with ``#CHR``, ``START`` (0-based), ``END``, and *id_col*.
    """
    genes = read_gtf(gtf_file, ("gene",))["gene"]
    if id_col != "gene_id":
        genes = genes.rename(columns={"gene_id": id_col})
    return genes


def write_sample_ids(
    sample_id: str,
    dataset_ids: list,
    sample_types: list,
    assay_types: list,
    out_file: str,
    rdr_base_dataset_ids=None,
):
    """Write ``sample_ids.tsv``, one row per observation of the bb matrices.

    Every column is a sample-file record key spelled the same way, except the derived
    leading ``SAMPLE``.

    Args:
        sample_id: Sample (patient) id, one per file.
        dataset_ids: Dataset id per observation, in matrix-observation order.
        sample_types: ``tumor``/``normal`` per observation.
        assay_types: Assay type per observation; the column is omitted when None.
        out_file: Output TSV path.
        rdr_base_dataset_ids: RDR baseline dataset id per observation; the column is
            omitted when None.

    Returns:
        The DataFrame written.
    """
    record = {
        "sample_id": sample_id,
        "dataset_id": dataset_ids,
        "sample_type": sample_types,
        "assay_type": assay_types,
        "rdr_base_dataset_id": rdr_base_dataset_ids,
    }
    # a multiome pair shares one dataset_id, so the assay is what keeps SAMPLE unique
    repeated = len(set(dataset_ids)) != len(dataset_ids)
    samples = [
        f"{sample_id}_{d}_{a}" if repeated else f"{sample_id}_{d}"
        for d, a in zip(dataset_ids, assay_types)
    ]
    assert len(set(samples)) == len(samples), (
        f"sample_ids.tsv, duplicate SAMPLE: {sorted(samples)}"
    )
    sample_dict = {"SAMPLE": samples}
    sample_dict.update({key: val for key, val in record.items() if val is not None})
    sample_df = pd.DataFrame(sample_dict)
    sample_df.to_csv(out_file, sep="\t", header=True, index=False)
    return sample_df


def write_bb_file(bbs: pd.DataFrame, out_file: str):
    """Write ``bb.tsv.gz``, the feature axis of every bb matrix.

    The three coordinate columns are required and the optional ones are written when
    present, so one schema covers all three modes; the binning internals (``bb_id``,
    ``BLOCKSIZE``, ``seg_id``, ``PS``, the cluster keys) are dropped.

    Args:
        bbs: bbs carrying at least ``#CHR``, ``START``, ``END``.
        out_file: Output TSV path; ``.gz`` is compressed by pandas.

    Returns:
        The DataFrame written.
    """
    bb_cols = ["#CHR", "START", "END"]
    missing = [c for c in bb_cols if c not in bbs.columns]
    assert not missing, f"bb table, missing column(s) {missing}"
    bb_cols += [
        c
        for c in ("#SNPS", "region_id", "switchprobs", "feature_id", "#feature")
        if c in bbs.columns
    ]
    bb_out = bbs[bb_cols]
    bb_out.to_csv(out_file, sep="\t", header=True, index=False)
    return bb_out


def write_snp_info(
    snps: pd.DataFrame,
    out_file: str,
):
    """Write the SNP feature axis of the allele matrices.

    ``PS`` (the phaser's phase-set label, which becomes the binning phase clusters) and
    ``seg_id`` are carried only when present.

    Args:
        snps: Filtered SNPs, in matrix-feature order.
        out_file: Output TSV path.

    Returns:
        The DataFrame written.
    """
    snp_cols = ["#CHR", "POS", "POS0", "START", "END", "GT", "PHASE"]
    # upstream phaser's phaseset label
    if "PS" in snps.columns:
        snp_cols.append("PS")
    logging.info(f"phase set (PS) column carried: {'PS' in snps.columns}")

    snp_cols += ["region_id"]
    if "seg_id" in snps.columns:
        snp_cols.append("seg_id")
    snp_cols += ["feature_id", "feature_type"]
    snp_info = snps[snp_cols]
    snp_info.to_csv(out_file, sep="\t", header=True, index=False)
    return snp_info
