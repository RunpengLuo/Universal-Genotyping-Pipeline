"""Tile segment.bed into fixed-width windows and annotate the bias covariates.

Last update: 2026-08-11

Inputs:
- aux_dir/segment.bed: tiled per row, so windows never cross segments
- reference: genome FASTA, for the GC column
- mappability_bed: optional BED, for the MAP column
- aux_dir/repliseq/{name}.{reference_version}.bedGraph: optional, for the REPLI column
- genome_size: chrom sizes TSV
Outputs:
- aux_dir/windows.bed.gz: fixed bins with ids and covariate columns
"""

import logging
import os

snakemake_handle = snakemake  # noqa: F821

from utils import (
    set_omp_threads,
    setup_logging,
    log_hist,
    maybe_path,
    sort_df_chr,
    strip_chr_prefix,
)

set_omp_threads(snakemake_handle)
setup_logging(snakemake_handle.log[0])

import numpy as np
import pandas as pd
from pybedtools import BedTool

from segmentation_utils import build_fixedwidth_bins
from io_utils import read_segment_bed, read_bedgraph
from range_utils import assign_range_to_range


inp = snakemake_handle.input
p = snakemake_handle.params
out_bed = snakemake_handle.output["window_bed"]


segment_bed = inp["segment_bed"]
genome_size = inp["genome_size"]
reference = inp["reference"]
chroms = list(p["chroms"])
input_nochr = p["input_nochr"]
logging.info(f"build_window_bed: window_size={p['window_size']}, {len(chroms)} chroms")

segments = read_segment_bed(segment_bed)
segments["#CHR"] = segments["#CHR"].astype(str)
segments[["START", "END"]] = segments[["START", "END"]].astype(np.int64)
segments = segments[segments["#CHR"].isin(chroms)].reset_index(drop=True)

windows = build_fixedwidth_bins(segments, int(p["window_size"]))
windows = sort_df_chr(windows, ch="#CHR", pos="START")
logging.info(f"tiled {len(windows)} windows on {windows['#CHR'].nunique()} chroms")

bed_windows = windows[["#CHR", "START", "END"]].copy()
if input_nochr:
    bed_windows["#CHR"] = bed_windows["#CHR"].map(strip_chr_prefix)

# GC: per-window GC fraction via pybedtools nucleotide_content
bt = BedTool.from_dataframe(bed_windows)
nuc = bt.nucleotide_content(fi=reference).to_dataframe(disable_auto_names=True)
assert len(nuc) == len(windows), (
    f"bedtools nuc, got {len(nuc)} rows, expected {len(windows)}"
)
windows["GC"] = nuc["5_pct_gc"].values
logging.info("annotated GC")

# MAP: per-window mean mappability via pybedtools map
mappability_bed = maybe_path(inp["mappability_bed"])
if mappability_bed:

    def _to_genome_style(feature):
        """Rename one pybedtools interval to the genome's naming."""
        nochr_chrom = strip_chr_prefix(feature.chrom)
        feature.chrom = nochr_chrom if input_nochr else f"chr{nochr_chrom}"
        return feature

    n_windows = len(windows)
    win_bed = bed_windows.copy()
    win_bed["_df_idx"] = np.arange(n_windows)
    wb = BedTool.from_dataframe(win_bed).sort(g=genome_size)
    # streams; bedtools resolves both operands against genome_size
    map_bt = (
        BedTool(mappability_bed).each(_to_genome_style).saveas().sort(g=genome_size)
    )
    map_cov = pd.read_csv(
        wb.map(b=map_bt, c=4, o="mean", g=genome_size).fn,
        sep="\t",
        header=None,
        names=["#CHR", "START", "END", "_df_idx", "MAP"],
        dtype={"#CHR": str, "START": int, "END": int, "_df_idx": int, "MAP": str},
    )
    assert len(map_cov) == n_windows, (
        f"bedtools map, got {len(map_cov)} rows, expected {n_windows}"
    )
    map_cov["MAP"] = (
        pd.to_numeric(map_cov["MAP"], errors="coerce").fillna(0.0).clip(0.0, 1.0)
    )
    windows["MAP"] = map_cov.sort_values("_df_idx")["MAP"].values
    logging.info(f"annotated MAP from {os.path.basename(mappability_bed)}")
else:
    logging.info("no mappability_bed; MAP skipped")

# REPLI: bin each Repli-seq bedGraph by midpoint, average across tracks
bedgraphs = list(inp["bedgraphs"])
if bedgraphs:
    n_win = len(windows)
    signal_sum = np.zeros(n_win, dtype=np.float64)
    signal_count = np.zeros(n_win, dtype=np.int64)
    win_ranges = windows[["#CHR", "START", "END"]].assign(_win=np.arange(n_win))
    for lf in bedgraphs:
        bg = read_bedgraph(lf, chroms=chroms)
        bg, _ = assign_range_to_range(
            bg, win_ranges, "_win", rule="midpoint", dropna=True
        )
        # bincount, not `arr[idx] +=`: several bedGraph rows can share one window
        hit = bg["_win"].to_numpy()
        signal_sum += np.bincount(hit, weights=bg["signal"].to_numpy(), minlength=n_win)
        signal_count += np.bincount(hit, minlength=n_win)
    with np.errstate(invalid="ignore", divide="ignore"):
        windows["REPLI"] = np.round(
            np.where(signal_count > 0, signal_sum / signal_count, np.nan), 6
        )
    logging.info(f"annotated REPLI from {len(bedgraphs)} bedGraph track(s)")
else:
    logging.info("no Repli-seq bedGraphs; REPLI skipped")

windows.to_csv(out_bed, sep="\t", index=False, compression="gzip")
logging.info(
    f"wrote {len(windows)} windows [{', '.join(windows.columns)}] -> {out_bed}"
)

log_hist(windows["END"] - windows["START"], "window length (bp)")
