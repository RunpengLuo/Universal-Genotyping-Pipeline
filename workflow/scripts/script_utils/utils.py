import gzip
import os
import sys
import logging

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../config"))
from const import *  # noqa: F401,F403


def setup_logging(log):
    """Configure root logger to write INFO-level messages with timestamps to a file.

    Parameters
    ----------
    log : str
        Path to the log file.
    """
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers.clear()
    fh = logging.FileHandler(log)
    fh.setLevel(logging.INFO)
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    root.addHandler(fh)


def maybe_path(x):
    """Return None if *x* is an empty list or None, otherwise return *x* unchanged.

    Useful for coercing Snakemake optional inputs to ``None``.
    """
    if x == [] or x is None:
        return None
    return x


def strip_chr_prefix(name):
    """Drop a leading ``chr`` from a chromosome name."""
    name = str(name)
    return name[3:] if name.lower().startswith("chr") else name


def add_chr_prefix(series):
    """Prepend ``chr`` to a chromosome Series that lacks it."""
    series = series.astype(str)
    if len(series) and series.iloc[0].lower().startswith("chr"):
        return series
    return "chr" + series


def require_matching_chr_style(bed_file, input_nochr, key):
    """Raise unless a BED names contigs the way genome_size does."""
    opener = gzip.open if str(bed_file).endswith(".gz") else open
    with opener(bed_file, "rt") as fh:
        first = next(
            (
                ln
                for ln in fh
                if ln.strip() and not ln.startswith(("#", "track", "browser"))
            ),
            "",
        )
    if not first:
        return
    bed_nochr = not first.split()[0].lower().startswith("chr")
    if bed_nochr != input_nochr:
        raise ValueError(
            f"{key} names contigs {'without' if bed_nochr else 'with'} a chr prefix "
            f"({first.split()[0]!r}), but genome_size uses the other style; "
            "bedtools resolves both against genome_size, so they must agree"
        )


def chrom_sort_key(chrom):
    """Genomic sort key: autosomes numerically (any count), then X, Y, M, then unknowns.

    Accepts ``chr``-prefixed or bare names, as str or int.
    """
    core = str(chrom)
    if core.lower().startswith("chr"):
        core = core[3:]
    if core.isdigit():
        return (0, int(core), "")
    special = {"X": 1, "Y": 2, "M": 3, "MT": 3}
    if core.upper() in special:
        return (1, special[core.upper()], "")
    return (2, 0, core)


def sort_chroms(chroms: list):
    """Sort chromosome names in genomic order. See :func:`chrom_sort_key`."""
    assert len(chroms) != 0
    return sorted((str(c) for c in chroms), key=chrom_sort_key)


def is_canonical_chrom(chrom):
    """True for autosomes and X/Y (the contigs kept in a window BED)."""
    core = str(chrom)
    if core.lower().startswith("chr"):
        core = core[3:]
    return core.isdigit() or core.upper() in ("X", "Y")


def sort_df_chr(df: pd.DataFrame, ch="#CHR", pos="POS"):
    """Sort a DataFrame by chromosome (genomic order) then by position, in-place.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with chromosome and position columns.
    ch : str
        Name of the chromosome column.
    pos : str
        Name of the position column.

    Returns
    -------
    pd.DataFrame
        The same DataFrame, sorted in-place.
    """
    # Cast to str so an int-typed #CHR (autosomes-only file) matches the string
    # categories below instead of silently becoming all-NaN.
    df[ch] = df[ch].astype(str)
    chs = sort_chroms(df[ch].unique().tolist())
    df[ch] = pd.Categorical(df[ch], categories=chs, ordered=True)
    df.sort_values(by=[ch, pos], inplace=True, ignore_index=True)
    return df
