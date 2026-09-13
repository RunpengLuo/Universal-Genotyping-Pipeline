"""Per-window capture-target coverage, the on/off-target split of a capture assay.

Last update: 2026-08-26

Inputs:
- aux_dir/windows.bed.gz (or the configured window_bed): the shared fixed bins
- target_bed: the capture kit's target intervals (BED3+), e.g. IDT xGen, Agilent SureSelect
Outputs:
- aux_dir/window.target.npz: `mat`, float32 on-target bp fraction in [0, 1], one row per
  window, aligned to `read_window_bed(window_bed, chroms=chroms)`

A WES library splits its reads into a captured fraction at tens to hundreds of x and an
off-target background near 0x. That split is a property of the window grid, not of any one
dataset, so it is computed once here and read by both `rd_correct` and `combine_counts`.
"""

import logging

snakemake_handle = snakemake

from utils import setup_logging

setup_logging(snakemake_handle.log[0])

import numpy as np

from io_utils import read_BED, read_window_bed
from range_utils import range_overlap_bp

# inputs
window_bed = snakemake_handle.input["window_bed"]
target_bed = snakemake_handle.input["target_bed"]

# parameters
chroms = list(snakemake_handle.params["chroms"])

# outputs
out_window_target = snakemake_handle.output["window_target"]

bin_df = read_window_bed(window_bed, chroms=chroms)
targets = read_BED(target_bed)
targets = targets[targets["#CHR"].isin(chroms)]
logging.info(
    f"{len(bin_df)} windows on {len(chroms)} chromosomes, "
    f"{len(targets)} target intervals from {target_bed}"
)

tgt_bp = range_overlap_bp(bin_df, targets)
window_bp = (bin_df["END"] - bin_df["START"]).to_numpy(dtype=np.int64)
frac = np.divide(
    tgt_bp, window_bp, where=window_bp > 0, out=np.zeros(len(bin_df), dtype=np.float64)
)

on = frac > 0
logging.info(
    f"on-target windows: {int(on.sum())}/{len(bin_df)} "
    f"({on.mean() * 100:.2f}%), {tgt_bp.sum() / 1e6:.2f} Mbp of target sequence"
)
logging.info(
    f"target fraction of an on-target window: min={frac[on].min():.3f} "
    f"median={np.median(frac[on]):.3f} mean={frac[on].mean():.3f}"
    if on.any()
    else "no window overlaps a target interval"
)
assert on.any(), (
    f"{target_bed}: no target interval overlaps any window on {chroms}; "
    "check the reference build and the chromosome naming"
)

np.savez_compressed(out_window_target, mat=frac.astype(np.float32))
logging.info(f"wrote per-window target fraction to {out_window_target}")
