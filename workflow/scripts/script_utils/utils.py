"""Shared helpers: thread limits, logging, chromosome naming and ordering.

Imported before numpy by every Snakemake ``script:`` entry point, so this module
must stay import-light: ``set_omp_threads`` only limits the BLAS/OpenMP runtimes
if it runs before they load. pandas is therefore imported inside the one function
that needs it.
"""

import os
import sys
import logging

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../config"))
from const import *  # noqa: F401,F403

_OMP_THREAD_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
)


def set_omp_threads(snakemake_handle):
    """Cap every BLAS/OpenMP runtime at the rule's thread count.

    Must be called before numpy (or anything importing it) is imported: the
    runtimes read these variables once, when their shared library loads.

    Args:
        snakemake_handle: The injected ``snakemake`` object.

    Returns:
        The thread count applied.
    """
    threads = int(getattr(snakemake_handle, "threads", 1))
    for var in _OMP_THREAD_VARS:
        os.environ[var] = str(threads)
    return threads


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


def logging_snakemake(msg):
    """Log to Snakemake's run log (``.snakemake/log``) instead of the console.

    Snakemake attaches two handlers to its ``snakemake.logging`` logger: a stream
    handler on stderr and a file handler on the run log. Handing the record to the file
    handler alone keeps workflow-parse diagnostics in the run log and off the terminal;
    ``logger.info`` would reach both, and ``print`` neither.

    Falls back to ``print`` when there is no file handler: outside Snakemake (e.g.
    ``resources/scripts/validate_sample_file.py``) and under ``--dryrun``, which writes
    no log file at all.

    Args:
        msg: Message text, emitted as one unformatted INFO record.
    """

    def logfile_handler():
        """Snakemake's run-log file handler, or None when there is no run."""
        try:
            from snakemake.logging import logger
        except ImportError:
            return None
        return next(
            (h for h in logger.handlers if isinstance(h, logging.FileHandler)), None
        )

    handler = logfile_handler()
    if handler is None:
        print(msg)
        return
    handler.handle(
        logging.LogRecord("snakemake", logging.INFO, __file__, 0, msg, None, None)
    )


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
    """Prepend ``chr`` to every chromosome of a Series that lacks it."""
    series = series.astype(str)
    return series.where(series.str.lower().str.startswith("chr"), "chr" + series)


def match_chr_style(name, input_nochr):
    """Rename one contig to the convention ``genome_size`` declares."""
    core = strip_chr_prefix(name)
    return core if input_nochr else f"chr{core}"


def chrom_sort_key(chrom):
    """Genomic sort key: autosomes numerically (any count), then X, Y, M, then unknowns.

    Accepts ``chr``-prefixed or bare names, as str or int.
    """
    core = strip_chr_prefix(chrom)
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
    core = strip_chr_prefix(chrom)
    return core.isdigit() or core.upper() in ("X", "Y")


def sort_df_chr(df, ch="#CHR", pos="POS"):
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
    import pandas as pd

    # Cast to str so an int-typed #CHR (autosomes-only file) matches the string
    # categories below instead of silently becoming all-NaN.
    df[ch] = df[ch].astype(str)
    chs = sort_chroms(df[ch].unique().tolist())
    df[ch] = pd.Categorical(df[ch], categories=chs, ordered=True)
    df.sort_values(by=[ch, pos], inplace=True, ignore_index=True)
    return df
