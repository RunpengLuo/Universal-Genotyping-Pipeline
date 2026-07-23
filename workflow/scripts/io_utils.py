import os
from collections import OrderedDict

import pandas as pd
import numpy as np

from utils import *


def get_chr_sizes(sz_file: str):
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


def get_standard_chroms(reference_version, genome_size_file):
    """Standard chromosomes (autosomes + X/Y) present in the genome-size file."""
    all_chroms = set(get_chr_sizes(genome_size_file).keys())
    candidates = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    if reference_version not in CHR_STYLE_REFVERS:
        candidates |= {str(c) for c in list(range(1, 23)) + ["X", "Y"]}
    return candidates & all_chroms


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
    snps["RAW_SNP_IDX"] = np.arange(len(snps))
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


def compute_depth_statistics(dp_raw, win_df, sample_ids):
    """Compute per-chromosome and whole-genome mean/median depth per sample.

    Returns a DataFrame with columns: SAMPLE, #CHR, mean_depth, median_depth.
    """
    chroms = win_df["#CHR"].to_numpy()
    sorted_chroms = sort_chroms(win_df["#CHR"].unique().tolist())
    rows = []
    for chrom in sorted_chroms:
        mask = chroms == chrom
        for s in range(len(sample_ids)):
            vals = dp_raw[mask, s]
            rows.append(
                [sample_ids[s], chrom, float(np.mean(vals)), float(np.median(vals))]
            )
    for s in range(len(sample_ids)):
        vals = dp_raw[:, s]
        rows.append(
            [sample_ids[s], "TOTAL", float(np.mean(vals)), float(np.median(vals))]
        )
    return pd.DataFrame(rows, columns=["SAMPLE", "#CHR", "mean_depth", "median_depth"])


def compute_snp_statistics(
    raw_snps_list,
    modalities,
    base_snps,
    chrom_sizes,
    chroms,
):
    """Compute per-chromosome SNP genotyping statistics.

    Returns a DataFrame with columns:
    #CHR, LENGTH, #SNPS-DNA, #SNPS-RNA, #SNPS-shared, #SNPS-union,
    #SNPS-het, #SNPS-hom_alt, #SNPS-hom_ref, #SNPS-dropped, #SNPS-kept.
    """
    chrom_list = [f"chr{c}" for c in chroms]
    modality_keys = {}
    for idx, modality in enumerate(modalities):
        modality_keys[modality] = set(raw_snps_list[idx]["KEY"])

    has_dna = "DNA" in modality_keys
    has_rna = "RNA" in modality_keys
    if has_dna and has_rna:
        shared_keys = modality_keys["DNA"] & modality_keys["RNA"]
    else:
        shared_keys = set()

    rows = []
    for chrom in chrom_list:
        length = chrom_sizes.get(chrom, 0)
        ch_mask = base_snps["#CHROM"] == chrom

        n_dna = 0
        n_rna = 0
        if has_dna:
            n_dna = int(
                (raw_snps_list[modalities.index("DNA")]["#CHROM"] == chrom).sum()
            )
        if has_rna:
            n_rna = int(
                (raw_snps_list[modalities.index("RNA")]["#CHROM"] == chrom).sum()
            )

        ch_keys = set(base_snps.loc[ch_mask, "KEY"])
        n_shared = len(ch_keys & shared_keys)
        n_union = int(ch_mask.sum())

        ch_snps = base_snps.loc[ch_mask]
        n_het = int((ch_snps["SAMPLE"] == "0/1").sum())
        n_hom_alt = int((ch_snps["SAMPLE"] == "1/1").sum())
        n_hom_ref = int((ch_snps["SAMPLE"] == "0/0").sum())
        n_dropped = n_union - n_het - n_hom_alt - n_hom_ref
        n_kept = n_het + n_hom_alt

        rows.append(
            [
                chrom,
                length,
                n_dna,
                n_rna,
                n_shared,
                n_union,
                n_het,
                n_hom_alt,
                n_hom_ref,
                n_dropped,
                n_kept,
            ]
        )

    # TOTAL row
    total_len = sum(chrom_sizes.get(c, 0) for c in chrom_list)
    total_dna = sum(r[2] for r in rows)
    total_rna = sum(r[3] for r in rows)
    total_shared = sum(r[4] for r in rows)
    total_union = sum(r[5] for r in rows)
    total_het = sum(r[6] for r in rows)
    total_hom_alt = sum(r[7] for r in rows)
    total_hom_ref = sum(r[8] for r in rows)
    total_dropped = sum(r[9] for r in rows)
    total_kept = sum(r[10] for r in rows)
    rows.append(
        [
            "TOTAL",
            total_len,
            total_dna,
            total_rna,
            total_shared,
            total_union,
            total_het,
            total_hom_alt,
            total_hom_ref,
            total_dropped,
            total_kept,
        ]
    )

    columns = [
        "#CHR",
        "LENGTH",
        "#SNPS-DNA",
        "#SNPS-RNA",
        "#SNPS-shared",
        "#SNPS-union",
        "#SNPS-het",
        "#SNPS-hom_alt",
        "#SNPS-hom_ref",
        "#SNPS-dropped",
        "#SNPS-kept",
    ]
    return pd.DataFrame(rows, columns=columns)


def _read_gtf(gtf_file: str, feature_type: str) -> pd.DataFrame:
    """Parse a GTF file and return records of the requested feature type.

    Returns a DataFrame with ``#CHR``, ``START`` (0-based), ``END``,
    and ``gene_id`` columns.
    """
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
        gtf["feature"] == feature_type, ["seqname", "start", "end", "attributes"]
    ]
    gene_ids = gtf["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
    return pd.DataFrame(
        {
            "#CHR": gtf["seqname"].values,
            "START": gtf["start"].values - 1,  # GTF is 1-based → 0-based
            "END": gtf["end"].values,  # GTF end is inclusive → half-open
            "gene_id": gene_ids.values,
        }
    )


def read_genes_gtf_file(gtf_file: str, id_col="gene_ids"):
    """Parse a GTF file and return gene-level records with genomic coordinates.

    Extracts gene features, deduplicates by ``gene_id``, and converts to
    0-based BED-like coordinates.

    Parameters
    ----------
    gtf_file : str
        Path to a GTF annotation file.
    id_col : str
        Column name for the gene identifier in the output.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START`` (0-based), ``END``, and *id_col*.
    """
    genes = _read_gtf(gtf_file, "gene")
    genes = genes.drop_duplicates("gene_id", keep="first")
    if id_col != "gene_id":
        genes = genes.rename(columns={"gene_id": id_col})
    return genes


def read_exons_gtf_file(gtf_file: str):
    """Parse a GTF file and return exon-level records with genomic coordinates.

    Extracts exon features and converts to 0-based BED-like coordinates.

    Parameters
    ----------
    gtf_file : str
        Path to a GTF annotation file.

    Returns
    -------
    pd.DataFrame
        DataFrame with ``#CHR``, ``START`` (0-based), ``END``, and ``gene_id``.
    """
    return _read_gtf(gtf_file, "exon")
