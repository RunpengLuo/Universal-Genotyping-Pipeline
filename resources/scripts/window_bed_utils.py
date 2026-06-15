"""Shared helpers for WGS and WES window BED generation.

Provides constants, chromosome handling, blacklist subtraction, region ID
assignment, and per-window covariate computation (GC, mappability, repli-seq).
"""

import glob
import os
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd
from pybedtools import BedTool

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

ALLOWED_REFVERS = ("hg19", "hg38", "chm13v2", "mm10")
CHR_STYLE_REFVERS = ("hg19", "hg38", "chm13v2", "mm10")
REPLISEQ_REFVERS = ("hg19", "hg38")

def chrom_sort_key(chrom):
    """Genomic sort key: autosomes numerically (any count), then X, Y."""
    core = str(chrom)
    if core.lower().startswith("chr"):
        core = core[3:]
    if core.isdigit():
        return (0, int(core))
    return (1, {"X": 1, "Y": 2}.get(core.upper(), 99))


def is_canonical_chrom(chrom):
    """True for autosomes and X/Y (the contigs kept in a window BED)."""
    core = str(chrom)
    if core.lower().startswith("chr"):
        core = core[3:]
    return core.isdigit() or core.upper() in ("X", "Y")

UCSC_REPLISEQ_BASE = (
    "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq"
)
REPLISEQ_BIGWIG_FILES = [
    "wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig",
    "wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig",
    "wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig",
]

LIFTOVER_CHAIN_URL = (
    "https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_chrom_sizes(genome_size_file):
    """Read chromosome sizes file into a DataFrame."""
    return pd.read_csv(
        genome_size_file,
        sep="\t",
        header=None,
        names=["chrom", "size"],
        dtype={"chrom": str, "size": np.int64},
    )


def _read_bed3(bed_file):
    """Read a BED file and return a DataFrame with #CHR, START, END."""
    return pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        comment="#",
        usecols=[0, 1, 2],
        names=["#CHR", "START", "END"],
        dtype={"#CHR": str, "START": np.int64, "END": np.int64},
    )


def _tile_region(chrom, start, end, window_size):
    """Tile a single region into fixed-size windows, merging last short bin.

    If the last window is shorter than window_size // 2 and there is a previous
    window, it is merged into the previous window (extending its END).  Regions
    smaller than window_size // 2 still emit one window.
    """
    rows = []
    pos = start
    while pos < end:
        w_end = min(pos + window_size, end)
        rows.append([chrom, pos, w_end])
        pos = w_end
    # Merge last bin into previous if undersized
    if len(rows) > 1 and (rows[-1][2] - rows[-1][1]) < window_size // 2:
        rows[-2][2] = rows[-1][2]
        rows.pop()
    return rows


def _sort_windows(windows):
    """Sort windows by genomic chromosome order then start; drop non-canonical contigs."""
    windows = windows[windows["#CHR"].map(is_canonical_chrom)].copy()
    windows["_chrom_rank"] = windows["#CHR"].map(chrom_sort_key)
    windows = windows.sort_values(["_chrom_rank", "START"]).reset_index(drop=True)
    return windows.drop(columns=["_chrom_rank"])


# ---------------------------------------------------------------------------
# Reference version helpers
# ---------------------------------------------------------------------------


def get_standard_chroms(reference_version, genome_size_file):
    """Derive valid standard chromosomes from reference version and genome size file."""
    all_chroms = set(_read_chrom_sizes(genome_size_file)["chrom"])
    candidates = {f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"]}
    if reference_version not in CHR_STYLE_REFVERS:
        candidates |= {str(c) for c in list(range(1, 23)) + ["X", "Y"]}
    return candidates & all_chroms


# ---------------------------------------------------------------------------
# Blacklist subtraction
# ---------------------------------------------------------------------------


def subtract_blacklist(windows, blacklist_bed_file):
    """Remove windows overlapping blacklisted regions."""
    bl = _read_bed3(blacklist_bed_file)
    keep = np.ones(len(windows), dtype=bool)
    for chrom in windows["#CHR"].unique():
        wm = (windows["#CHR"] == chrom).to_numpy()
        rm = (bl["#CHR"] == chrom).to_numpy()
        if not rm.any():
            continue
        w_starts = windows.loc[wm, "START"].to_numpy()
        w_ends = windows.loc[wm, "END"].to_numpy()
        bl_starts = bl.loc[rm, "START"].to_numpy()
        bl_ends = bl.loc[rm, "END"].to_numpy()
        sort_idx = np.argsort(bl_starts)
        bl_starts, bl_ends = bl_starts[sort_idx], bl_ends[sort_idx]
        right_bounds = np.searchsorted(bl_starts, w_ends, side="left")
        overlaps = np.array(
            [
                np.any(bl_ends[: right_bounds[i]] > w_starts[i])
                for i in range(len(w_starts))
            ]
        )
        keep[np.where(wm)[0]] = ~overlaps
    return windows.loc[keep].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Region ID assignment
# ---------------------------------------------------------------------------


def assign_region_id(windows, region_bed_file):
    """Assign region_id to each window by mapping midpoints to region intervals.

    Uses the BED's 4th column as ``region_id`` when present; otherwise falls
    back to ``CHR:START-END``.
    """
    regions = pd.read_csv(
        region_bed_file, sep="\t", header=None, comment="#", dtype={0: str}
    )
    assert len(regions.columns) >= 3, "invalid regions BED format"
    columns = ["#CHR", "START", "END"]
    if len(regions.columns) > 3:
        regions = regions.iloc[:, :4].copy()
        columns += ["region_id"]
    regions.columns = columns
    regions[["START", "END"]] = regions[["START", "END"]].astype(np.int64)
    if "region_id" not in regions.columns:
        regions["region_id"] = (
            regions["#CHR"] + ":" + regions["START"].astype(str)
            + "-" + regions["END"].astype(str)
        )

    mids = (windows["START"] + windows["END"]) // 2
    result = pd.Series(pd.NA, index=windows.index, dtype="object")

    for chrom in regions["#CHR"].unique():
        win_mask = windows["#CHR"] == chrom
        if not win_mask.any():
            continue
        reg_ch = regions.loc[regions["#CHR"] == chrom].sort_values("START")
        starts = reg_ch["START"].to_numpy()
        ends = reg_ch["END"].to_numpy()
        ids = reg_ch["region_id"].to_numpy()
        positions = mids[win_mask].to_numpy()

        idx = np.searchsorted(starts, positions, side="right") - 1
        safe_idx = idx.clip(min=0)
        valid = (idx >= 0) & (positions < ends[safe_idx])
        win_indices = windows.index[win_mask]
        result.loc[win_indices[valid]] = ids[idx[valid]]

    return result


# ---------------------------------------------------------------------------
# Per-window covariates
# ---------------------------------------------------------------------------


def compute_gc(windows_df, reference):
    """Compute per-window GC fraction via pybedtools nucleotide_content."""
    bt = BedTool.from_dataframe(windows_df[["#CHR", "START", "END"]])
    nuc = bt.nucleotide_content(fi=reference).to_dataframe(disable_auto_names=True)
    windows_df = windows_df.copy()
    windows_df["GC"] = nuc["5_pct_gc"].values
    return windows_df


def compute_mappability(windows_df, mappability_file, genome_size_file):
    """Compute per-window mean mappability via pybedtools map."""
    n_windows = len(windows_df)
    win_bed = windows_df[["#CHR", "START", "END"]].copy()
    win_bed["_idx"] = np.arange(n_windows)
    bt = BedTool.from_dataframe(win_bed).sort(g=genome_size_file)
    map_bt = BedTool(mappability_file).sort(g=genome_size_file)
    result_bt = bt.map(b=map_bt, c=4, o="mean", g=genome_size_file)
    # Read raw TSV instead of to_dataframe() which can silently drop rows
    map_cov = pd.read_csv(
        result_bt.fn,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "_idx", "MAP"],
        dtype={"#CHR": str, "START": int, "END": int, "_idx": int, "MAP": str},
    )
    map_cov["MAP"] = (
        pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    )
    assert len(map_cov) == n_windows, (
        f"bedtools map returned {len(map_cov)} rows, expected {n_windows}"
    )
    map_cov = map_cov.sort_values("_idx")
    return map_cov["MAP"].values


# ---------------------------------------------------------------------------
# Replication timing (Repli-seq)
# ---------------------------------------------------------------------------


def _download(url, dest):
    """Download a file with wget."""
    subprocess.check_call(["wget", "-q", "-O", dest, url])


def _bigwig_to_bedgraph(bigwig, bedgraph):
    """Convert bigWig to bedGraph using UCSC bigWigToBedGraph."""
    subprocess.check_call(["bigWigToBedGraph", bigwig, bedgraph])


def _liftover(input_bed, chain, output_bed, unmapped):
    """Run UCSC liftOver."""
    subprocess.check_call(["liftOver", input_bed, chain, output_bed, unmapped])


def _prepare_lifted_bedgraphs(bigwig_dir, work_dir, chain):
    """Download bigWigs if needed, convert to bedGraph, liftOver to hg38.

    Returns list of lifted bedGraph file paths.
    """
    os.makedirs(bigwig_dir, exist_ok=True)
    n_files = len(REPLISEQ_BIGWIG_FILES)
    for i, fname in enumerate(REPLISEQ_BIGWIG_FILES, 1):
        dest = os.path.join(bigwig_dir, fname)
        if not os.path.isfile(dest):
            print(f"\r  downloading bigWigs [{i}/{n_files}]", end="", flush=True)
            _download(f"{UCSC_REPLISEQ_BASE}/{fname}", dest)
    print()

    bigwig_files = sorted(glob.glob(os.path.join(bigwig_dir, "*WaveSignal*.bigWig")))
    if not bigwig_files:
        sys.exit(f"Error: no WaveSignal bigWig files found in {bigwig_dir}")

    lifted_dir = os.path.join(work_dir, "lifted")
    os.makedirs(lifted_dir, exist_ok=True)

    lifted_files = []
    for i, bw in enumerate(bigwig_files, 1):
        name = os.path.basename(bw).replace(".bigWig", "")
        bg_file = os.path.join(lifted_dir, f"{name}.hg19.bedGraph")
        lifted_file = os.path.join(lifted_dir, f"{name}.hg38.bedGraph")
        unmapped_file = os.path.join(lifted_dir, f"{name}.unmapped")

        if not os.path.isfile(lifted_file):
            print(f"\r  liftOver [{i}/{len(bigwig_files)}]", end="", flush=True)
            _bigwig_to_bedgraph(bw, bg_file)
            _liftover(bg_file, chain, lifted_file, unmapped_file)
            os.remove(bg_file)
        lifted_files.append(lifted_file)
    print()
    return lifted_files


def _bin_bedgraph_signals(windows_df, lifted_files, standard_chroms):
    """Bin lifted bedGraph signals into windows. Returns per-window mean signal."""
    n_win = len(windows_df)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int32)

    for i, lf in enumerate(lifted_files, 1):
        print(f"\r  binning repli-seq [{i}/{len(lifted_files)}]", end="", flush=True)
        bg = pd.read_csv(
            lf,
            sep="\t",
            header=None,
            names=["chrom", "start", "end", "signal"],
            dtype={
                "chrom": str,
                "start": np.int64,
                "end": np.int64,
                "signal": np.float64,
            },
        )
        bg = bg[bg["chrom"].isin(standard_chroms)].reset_index(drop=True)

        for chrom, grp in bg.groupby("chrom", sort=False):
            win_mask = windows_df["#CHR"] == chrom
            win_idx = np.where(win_mask)[0]
            if len(win_idx) == 0:
                continue

            win_starts = windows_df["START"].to_numpy()[win_idx]
            sort_order = np.argsort(win_starts)
            win_starts = win_starts[sort_order]
            win_idx = win_idx[sort_order]

            bg_mids = ((grp["start"] + grp["end"]) // 2).to_numpy()
            bin_idx = np.searchsorted(win_starts, bg_mids, side="right") - 1
            valid = (bin_idx >= 0) & (bin_idx < len(win_idx))
            global_idx = win_idx[bin_idx[valid]]
            signal_sum[global_idx] += grp["signal"].to_numpy()[valid]
            signal_count[global_idx] += 1
    print()

    with np.errstate(invalid="ignore", divide="ignore"):
        mean_signal = np.where(signal_count > 0, signal_sum / signal_count, np.nan)

    n_covered = int((signal_count > 0).sum())
    print(f"  {n_covered}/{n_win} windows have replication timing data")
    return np.round(mean_signal, 6)


def compute_repliseq(windows_df, standard_chroms, args):
    """Compute per-window mean replication timing from ENCODE Repli-seq bigWigs.

    Handles work_dir/chain setup, then delegates to helpers.
    """
    if args.work_dir:
        work_dir = args.work_dir
        os.makedirs(work_dir, exist_ok=True)
        cleanup = False
    else:
        work_dir = tempfile.mkdtemp(prefix="repliseq_")
        cleanup = True

    chain = args.chain
    if chain is None:
        chain = os.path.join(work_dir, "hg19ToHg38.over.chain.gz")
        if not os.path.isfile(chain):
            print("  Downloading chain file ...")
            _download(LIFTOVER_CHAIN_URL, chain)
    if not os.path.isfile(chain):
        sys.exit(f"Error: chain file not found: {chain}")

    bigwig_dir = args.bigwig_dir or os.path.join(work_dir, "bigwigs")
    lifted_files = _prepare_lifted_bedgraphs(bigwig_dir, work_dir, chain)
    result = _bin_bedgraph_signals(windows_df, lifted_files, standard_chroms)

    if cleanup:
        import shutil

        shutil.rmtree(work_dir)

    return result


# ---------------------------------------------------------------------------
# Summary reporting
# ---------------------------------------------------------------------------


def print_summary(windows):
    """Print per-window and per-chromosome summary statistics."""
    sizes = windows["END"] - windows["START"]

    print()
    print("  Summary:")
    print(f"    {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")
    print(f"    {windows['region_id'].nunique()} distinct regions")
    print(
        f"    window size: min={sizes.min()}, median={int(sizes.median())}, "
        f"max={sizes.max()}, mean={sizes.mean():.0f}"
    )
    total_bp = int(sizes.sum())
    print(f"    total coverage: {total_bp:,} bp ({total_bp / 1e9:.2f} Gb)")

    # Per-chromosome table
    header = f"    {'chrom':<8} {'windows':>8} {'regions':>8} {'coverage_Mb':>12} {'median_size':>12}"
    print(header)

    present_chroms = sorted(
        (c for c in windows["#CHR"].unique() if is_canonical_chrom(c)),
        key=chrom_sort_key,
    )
    for chrom in present_chroms:
        mask = windows["#CHR"] == chrom
        if not mask.any():
            continue
        ch_sizes = sizes[mask]
        ch_regions = windows.loc[mask, "region_id"].nunique()
        row = (
            f"    {chrom:<8} {int(mask.sum()):>8} {ch_regions:>8} "
            f"{ch_sizes.sum() / 1e6:>12.2f} {int(ch_sizes.median()):>12}"
        )
        print(row)
    print()
