import logging
import os
import tempfile
from collections import OrderedDict

import pandas as pd
import numpy as np

from const import (
    GTF_COLUMNS,
    RANGER_MATRIX_H5,
    RANGER_SPATIAL_DIR,
    SAMPLE_ID_COLNAMES,
)
from utils import add_chr_prefix, sort_chroms, sort_df_chr


def read_chrom_sizes(sz_file: str):
    """Read a two-column chromosome-sizes file and return an OrderedDict mapping name to length.

    Parameters
    ----------
    sz_file : str
        Path to a tab-separated file with columns (chromosome, size).

    Returns
    -------
    OrderedDict[str, int]
        Chromosome name to integer length.
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
    """
    load VCF file as dataframe.
    If phased, parse GT[0] as USEREF, check PS
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

    Input is the output of ``pileup_snps_bulk_bcftools``: tab-separated
    ``#CHROM POS REF ALT AD``, where ALT is a comma-list (e.g. ``C,<*>``) and AD is a
    comma-list of allele depths (``ref,alt1,...``), one row per het locus with reads.

    Only the PARENT ALT allele is counted: its depth is looked up in the locus's own ALT
    list and is 0 when that allele was not observed. ``DP = ref + alt`` is therefore the
    usable het depth, so a downstream ``REF = DP - ALT`` is exact.

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


def read_snp_mats_bulk(snp_info_file, tot_file, a_file, b_file):
    """Read the joint bulk SNP table and dense T/A/B matrices, genomically sorted.

    Returns ``(snps, tot_mtx, a_mtx, b_mtx)`` with SNP rows in ``#CHR``/``POS0``
    order and the matrices permuted to match.
    """
    snps = pd.read_table(snp_info_file, sep="\t")
    tot_mtx = np.load(tot_file)["mat"].astype(np.int32)
    a_mtx = np.load(a_file)["mat"].astype(np.int32)
    b_mtx = np.load(b_file)["mat"].astype(np.int32)

    snps["_row"] = np.arange(len(snps))
    snps = sort_df_chr(snps, ch="#CHR", pos="POS0").reset_index(drop=True)
    perm = snps["_row"].to_numpy()
    tot_mtx, a_mtx, b_mtx = tot_mtx[perm], a_mtx[perm], b_mtx[perm]
    snps = snps.drop(columns="_row")
    return snps, tot_mtx, a_mtx, b_mtx


def read_BED(bed_file: str, addchr=True, extra_columns=("region_id", "seg_id")):
    """Read a BED file: the first 3 columns are ``#CHR``/``START``/``END``.

    Any further columns are named from ``extra_columns`` in order. When an
    ``extra_columns`` entry has no column in the file it is filled: ``region_id``
    falls back to ``CHR:START-END`` and ``seg_id`` to ``region_id`` (the
    build_segment_bed default). Pass ``extra_columns=()`` for a plain BED3
    (e.g. a blacklist).
    """
    df = pd.read_table(bed_file, sep="\t", header=None, dtype={0: "string"})
    assert len(df.columns) >= 3, (
        f"BED file, expected >=3 columns, got {len(df.columns)}"
    )
    n_extra = min(len(df.columns) - 3, len(extra_columns))
    columns = ["Chromosome", "Start", "End"] + list(extra_columns[:n_extra])
    df = df.iloc[:, : 3 + n_extra].copy()
    df.columns = columns
    if not str(df["Chromosome"].iloc[0]).startswith("chr") and addchr:
        df["Chromosome"] = "chr" + df["Chromosome"].astype(str)
    df["#CHR"] = df["Chromosome"]
    df["START"] = df["Start"]
    df["END"] = df["End"]

    if "region_id" in extra_columns and "region_id" not in df.columns:
        df["region_id"] = (
            df["#CHR"] + ":" + df["START"].astype(str) + "-" + df["END"].astype(str)
        )
    if "seg_id" in extra_columns and "seg_id" not in df.columns:
        df["seg_id"] = df["region_id"]
    return df


def read_bedgraph(bg_file: str, chroms=None):
    """Read a bedGraph track: ``chrom start end value``, 0-based half-open.

    Contig names are chr-normalized on ingest, like every other reader here, so a
    track that spells its contigs bare still matches *chroms*.

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
    """Read a barcode file (one barcode per line) and return as a list of strings.

    Parameters
    ----------
    bc_file : str
        Path to a text file with one barcode per line.

    Returns
    -------
    list[str]
        Barcodes.
    """
    barcodes = (
        pd.read_table(bc_file, sep="\t", header=None, dtype=str).iloc[:, 0].tolist()
    )
    return barcodes


def read_full_barcodes(path: str):
    """Read a 2-column REP_ID,BARCODE TSV (with header) into a DataFrame."""
    return pd.read_table(path, sep="\t", header=0, dtype=str)


def read_chunks_from_atac_fragments(frag_file: str, chunksize=5_000_000):
    """Read a 10x ATAC fragment file in chunks.

    Each record of ``atac_fragments.tsv.gz`` is one deduplicated fragment, columns
    ``chrom, chromStart, chromEnd, barcode, readSupport, strand``; only the first four
    are read. One sample runs to hundreds of millions of records, hence the chunking.
    Contig names are left as the file spells them, so a caller that filters rows first
    can chr-normalize the subset rather than every record.

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

    squidpy takes a directory, while the sample file names each spatial file
    individually so remote files can be fetched. The Space Ranger layout it expects is
    rebuilt as symlinks in a temporary directory - the feature matrix at the root, the
    rest under spatial/ - which lives only for the read. Names come from RANGER_* in
    const.py.

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

    Contig names are normalized to chr-notation; coordinates become 0-based
    half-open. Genes are deduplicated by ``gene_id``; other feature types keep
    every record.

    Args:
        gtf_file: Path to a GTF annotation file (optionally gzipped).
        feature_types: Feature types to extract, e.g. ``("gene", "exon")``.

    Returns:
        ``{feature_type: DataFrame}`` with ``#CHR``, ``START``, ``END``, ``gene_id``.
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
    colnames=None,
):
    """Write ``sample_ids.tsv``, one row per observation of the bb matrices.

    Every column but the leading ``SAMPLE`` is a sample-file record key renamed through
    ``SAMPLE_ID_COLNAMES``; ``SAMPLE`` is derived (``{sample_id}_{dataset_id}``) and has
    no record key. Column order follows the argument order.

    Args:
        sample_id: Sample (patient) id, one per file.
        dataset_ids: Dataset id per observation, in matrix-observation order.
        sample_types: ``tumor``/``normal`` per observation.
        assay_types: Assay type per observation; the column is omitted when None.
        out_file: Output TSV path.
        rdr_base_dataset_ids: RDR baseline dataset id per observation; the column is
            omitted when None.
        colnames: Record key -> column name overrides on top of ``SAMPLE_ID_COLNAMES``.

    Returns:
        The DataFrame written.
    """
    cols = {**SAMPLE_ID_COLNAMES, **(colnames or {})}
    record = {
        "sample_id": sample_id,
        "dataset_id": dataset_ids,
        "sample_type": sample_types,
        "assay_type": assay_types,
        "rdr_base_dataset_id": rdr_base_dataset_ids,
    }
    sample_dict = {"SAMPLE": [f"{sample_id}_{d}" for d in dataset_ids]}
    sample_dict.update(
        {cols[key]: val for key, val in record.items() if val is not None}
    )
    sample_df = pd.DataFrame(sample_dict)
    sample_df.to_csv(out_file, sep="\t", header=True, index=False)
    return sample_df
