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
from const import is_url

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

    Falls back to ``print`` when there is no file handler: outside Snakemake, and under
    ``--dryrun``, which writes no log file at all.

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


def log_hist(values, label, bins=20, width=48, fmt=".4g"):
    """Log a one-line summary plus a plain-ASCII histogram of *values*.

    For a distribution worth eyeballing but not worth a PDF. Bars are scaled to the
    tallest bin, so the shape is readable but the heights are relative; the counts are
    printed alongside. numpy is imported here, not at module scope, because this module
    loads before ``set_omp_threads`` has capped the BLAS thread count.

    All-integer values (counts) get integer bin edges, one bucket per value while the
    range fits in *bins*, so a count histogram has no fractional edge and no empty
    bucket between two attainable values.

    Args:
        values: 1-D numeric sequence; non-finite entries are dropped.
        label: Name of the quantity, used in the summary line.
        bins: Number of histogram bins, or the most for integer values.
        width: Character width of the tallest bar.
        fmt: Format spec for the summary statistics and bin edges.
    """
    import numpy as np

    v = np.asarray(values, dtype=float).ravel()
    v = v[np.isfinite(v)]
    if v.size == 0:
        logging.info(f"{label}: no finite values")
        return
    logging.info(
        f"{label}: n={v.size}  min={v.min():{fmt}}  median={np.median(v):{fmt}}  "
        f"mean={v.mean():{fmt}}  max={v.max():{fmt}}"
    )
    lo, hi = v.min(), v.max()
    step = 0
    if hi == lo:
        bins = 1
    elif np.all(v == np.floor(v)):
        step = max(1, int(np.ceil((hi - lo + 1) / bins)))
        bins = np.arange(lo, hi + step + 1, step)
    counts, edges = np.histogram(v, bins=bins)
    peak = counts.max()
    for i, count in enumerate(counts):
        bar = "#" * round(width * count / peak) if peak else ""
        if step == 1:
            logging.info(f"  {edges[i]:>22{fmt}} {count:>9d} |{bar}")
            continue
        # np.histogram closes only the last bin on the right
        close = "]" if i == len(counts) - 1 else ")"
        logging.info(
            f"  [{edges[i]:>10{fmt}}, {edges[i + 1]:>10{fmt}}{close} {count:>9d} |{bar}"
        )


def maybe_path(x):
    """Return None if *x* is an empty list or None, otherwise return *x* unchanged.

    Useful for coercing Snakemake optional inputs to ``None``.
    """
    if x == [] or x is None:
        return None
    return x


def check_local_path(path, label):
    """Assert an input path is usable: a URL by scheme, a local path by existence.

    A URL is never contacted. ``storage()`` fetches it when a job needs it, and a
    reachability probe here would put a network request in every dry run.

    Args:
        path: Local path or URL.
        label: Names where the path came from, e.g. ``"d1: files.alignment"``.

    Returns:
        *path* unchanged.

    Raises:
        AssertionError: *path* is local and does not exist.
    """
    if is_url(path):
        return path
    assert os.path.exists(path), f"{label}, path does not exist: {path}"
    return path


def strip_chr_prefix(name):
    """Drop a leading ``chr`` from a chromosome name."""
    name = str(name)
    return name[3:] if name.lower().startswith("chr") else name


def add_chr_prefix(series):
    """Prepend ``chr`` to every chromosome of a Series that lacks it."""
    series = series.astype(str)
    return series.where(series.str.lower().str.startswith("chr"), "chr" + series)


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
    assert len(chroms) != 0, "chromosome list is empty"
    return sorted((str(c) for c in chroms), key=chrom_sort_key)


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
