from collections import OrderedDict

import pandas as pd
import numpy as np

from const import GTF_COLUMNS
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
    assert ncols == 8 or ncols >= 10, "invalid VCF file"
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


def read_bcftools_counts(tsv_file: str):
    """Read a bcftools per-locus AD table into a DataFrame.

    Input is the output of ``pileup_snps_bulk_bcftools``: tab-separated
    ``#CHROM POS REF ALT AD`` where ALT is a comma-list (e.g. ``C,<*>``) and AD is a
    comma-list of allele depths (``ref,alt1,...``), one row per het locus with reads.

    Args:
        tsv_file: Path to the (optionally gzipped) counts TSV.

    Returns:
        DataFrame with columns ``#CHROM``, ``POS``, ``REF``, ``ALT`` (list[str]),
        ``AD`` (list[int]), ``KEY`` (``#CHROM_POS``, matching ``read_VCF``), and
        ``RAW_SNP_DF_IDX`` (file row order). Empty DataFrame if the file has no records.
    """
    df = pd.read_csv(
        tsv_file,
        sep="\t",
        header=None,
        names=["#CHROM", "POS", "REF", "ALT", "AD"],
        dtype={"#CHROM": "string", "REF": "string", "ALT": "string", "AD": "string"},
    )
    if df.empty:
        return df
    df["POS"] = df["POS"].astype(np.int64)
    if not str(df["#CHROM"].iloc[0]).startswith("chr"):
        df["#CHROM"] = "chr" + df["#CHROM"].astype(str)
    df["#CHROM"] = df["#CHROM"].str.replace("^chrMT$", "chrM", regex=True)
    df["KEY"] = df["#CHROM"].astype(str) + "_" + df["POS"].astype(str)
    df["RAW_SNP_DF_IDX"] = np.arange(len(df))
    df["ALT"] = df["ALT"].str.split(",")
    df["AD"] = df["AD"].str.split(",").apply(lambda xs: [int(x) for x in xs])
    return df


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
    assert len(df.columns) >= 3, "invalid BED format"
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


def read_BEDPE(bedpe_file: str, addchr=True):
    """Read a BEDPE of SV junctions (bedtools 10-column layout, 0-based like BED).

    Columns: ``#CHR1 START1 END1 #CHR2 START2 END2 [NAME SCORE STRAND1 STRAND2]``;
    only the first 6 are required. Both chrom columns are chr-normalized so they
    match a chr-prefixed region BED. An empty file returns an empty DataFrame.
    """
    names = [
        "#CHR1",
        "START1",
        "END1",
        "#CHR2",
        "START2",
        "END2",
        "NAME",
        "SCORE",
        "STRAND1",
        "STRAND2",
    ]
    try:
        df = pd.read_csv(
            bedpe_file, sep="\t", header=None, comment="#", dtype={0: str, 3: str}
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=names[:6])
    n = min(df.shape[1], len(names))
    df = df.iloc[:, :n].copy()
    df.columns = names[:n]
    for c in ("#CHR1", "#CHR2"):
        s = df[c].astype(str)
        df[c] = s if not addchr else s.where(s.str.startswith("chr"), "chr" + s)
    df["START1"] = df["START1"].astype(np.int64)
    df["START2"] = df["START2"].astype(np.int64)
    return df


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
