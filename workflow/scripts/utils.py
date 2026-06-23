import os
import sys
import logging

import pandas as pd
import numpy as np

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


def sort_chroms(chromosomes: list):
    """Sort chromosome names in genomic order. See :func:`chrom_sort_key`."""
    assert len(chromosomes) != 0
    return sorted((str(c) for c in chromosomes), key=chrom_sort_key)


def adaptive_dot_size(n_points, s_base=4, s_min=0.5, s_max=10, n_ref=5000):
    """Scale dot size inversely with point count.

    At *n_ref* points the size equals *s_base*; fewer points -> bigger dots,
    more points -> smaller dots, clamped to [s_min, s_max].
    """
    if n_points <= 0:
        return s_base
    return float(np.clip(s_base * n_ref / n_points, s_min, s_max))



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
