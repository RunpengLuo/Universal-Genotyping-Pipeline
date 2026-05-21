#!/usr/bin/env python3
"""Build a per-window bias correction BED for WGS data.

Tiles fixed-size windows within accessible regions, then computes per-window
GC content with optional mappability and replication timing columns.

Pipeline:
  1. Subtract --blacklist_bed from --region_bed (if provided)
  2. Tile directly within each remaining interval (every window fully contained)
  3. Merge undersized trailing windows into the previous window
  4. Assign region_id, compute GC/MAP/REPLI
  5. Sort and write output

Output is a gzipped TSV:
  #CHR  START  END  region_id  GC  [MAP]  [REPLI]

Dependencies:
  - pybedtools, pandas, numpy
  - bigWigToBedGraph, liftOver (UCSC tools) -- only if --repliseq
"""

import argparse
import os
import sys
import time

import pandas as pd
from pybedtools import BedTool

from window_bed_utils import (
    REPLISEQ_REFVERS,
    _sort_windows,
    _tile_region,
    assign_region_id,
    compute_gc,
    compute_mappability,
    compute_repliseq,
    get_standard_chroms,
    print_summary,
)


# ---------------------------------------------------------------------------
# Window generation
# ---------------------------------------------------------------------------


def generate_wgs_windows(
    genome_size_file, window_size, standard_chroms, region_bed, blacklist_bed
):
    """Tile windows directly within accessible regions.

    Reads the region BED, subtracts the blacklist, then tiles each remaining
    interval into fixed-size windows (merging undersized trailing windows).
    Every output window is fully contained in a blacklist-subtracted region.
    """
    regions_bt = BedTool(region_bed)
    if blacklist_bed:
        regions_bt = regions_bt.subtract(BedTool(blacklist_bed))

    reg_df = pd.read_csv(
        regions_bt.fn,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"],
        dtype={"Chromosome": str},
    )
    reg_df = reg_df[reg_df["Chromosome"].isin(standard_chroms)].reset_index(drop=True)

    rows = []
    for _, r in reg_df.iterrows():
        rows.extend(_tile_region(r["Chromosome"], r["Start"], r["End"], window_size))

    windows = pd.DataFrame(rows, columns=["#CHR", "START", "END"])
    sizes = windows["END"] - windows["START"]
    print(f"  {len(windows)} windows across {windows['#CHR'].nunique()} chromosomes")
    print(
        f"  window size: min={sizes.min()}, mean={sizes.mean():.0f}, "
        f"median={int(sizes.median())}, max={sizes.max()}"
    )
    return windows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build a per-window bias correction BED for WGS. "
            "Tiles fixed-size windows within accessible regions, then computes "
            "per-window GC content with optional mappability and replication timing."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Example:\n"
            "  python build_wgs_window_bed.py \\\n"
            "    --reference_version hg38 \\\n"
            "    --reference /path/to/hg38.fa \\\n"
            "    --genome_size /path/to/hg38.genome \\\n"
            "    --region_bed /path/to/regions.bed \\\n"
            "    --out_file gc.1kbp.hg38.bed.gz"
        ),
    )
    # Required
    parser.add_argument(
        "--reference_version",
        required=True,
        help="Reference genome version (hg19, hg38, chm13v2, mm10).",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference FASTA (must have .fai index).",
    )
    parser.add_argument(
        "--genome_size",
        required=True,
        help="Path to chromosome sizes file (tab-separated: chrom, size).",
    )
    parser.add_argument(
        "--region_bed",
        required=True,
        help="Accessible regions BED -- keep only windows overlapping these regions.",
    )
    parser.add_argument(
        "--out_file",
        required=True,
        help="Output gzipped TSV path (e.g. windows.1kbp.hg38.bed.gz).",
    )
    # Optional filtering
    parser.add_argument(
        "--blacklist_bed",
        default=None,
        help="Optional blacklist BED -- subtract these regions.",
    )
    # WGS options
    parser.add_argument(
        "--window", type=int, default=1000, help="Window size in bp (default: 1000)."
    )
    # Covariate options
    parser.add_argument(
        "--mappability_bed",
        default=None,
        help="Optional BED-format mappability track (4th column = score).",
    )
    parser.add_argument(
        "--repliseq",
        action="store_true",
        help="Enable replication timing from ENCODE Repli-seq bigWigs.",
    )
    parser.add_argument(
        "--chain",
        default=None,
        help="hg19-to-hg38 liftOver chain file (downloaded if omitted with --repliseq).",
    )
    parser.add_argument(
        "--bigwig_dir",
        default=None,
        help="Directory for pre-downloaded WaveSignal bigWig files.",
    )
    parser.add_argument(
        "--work_dir",
        default=None,
        help="Working directory for intermediate files (default: temp dir).",
    )
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def main():
    args = parse_args()

    # --- Validate inputs ---
    for name, path in [
        ("reference", args.reference),
        ("genome_size", args.genome_size),
        ("region_bed", args.region_bed),
        ("blacklist_bed", args.blacklist_bed),
        ("mappability_bed", args.mappability_bed),
    ]:
        if path and not os.path.isfile(path):
            sys.exit(f"Error: {name} not found: {path}")

    t0 = time.time()

    print("=== build_wgs_window_bed.py ===")
    print(f"  Ref version: {args.reference_version}")
    print(f"  Reference:   {args.reference}")
    print(f"  Genome size: {args.genome_size}")
    print(f"  Region BED:  {args.region_bed}")
    print(f"  Blacklist:   {args.blacklist_bed or 'none'}")
    print(f"  Window:      {args.window} bp")
    print(f"  Mappability: {args.mappability_bed or 'none'}")
    print(f"  Repli-seq:   {'yes' if args.repliseq else 'no'}")
    print(f"  Output:      {args.out_file}")
    print()

    standard_chroms = get_standard_chroms(args.reference_version, args.genome_size)
    print(f"  {len(standard_chroms)} standard chromosomes")

    # --- Step 1: Generate windows ---
    print("[1/6] WGS: tiling windows within regions ...")
    windows = generate_wgs_windows(
        args.genome_size,
        args.window,
        standard_chroms,
        args.region_bed,
        args.blacklist_bed,
    )

    # --- Step 2: Assign region_id, drop windows outside regions ---
    print("[2/6] Assigning region_id ...")
    windows["region_id"] = assign_region_id(windows, args.region_bed)
    n_before = len(windows)
    windows = windows[windows["region_id"].notna()].reset_index(drop=True)
    n_dropped = n_before - len(windows)
    assert n_dropped == 0, (
        f"{n_dropped}/{n_before} windows could not be assigned a region_id"
    )
    print(f"  {len(windows)} assigned")

    # --- Step 3: Compute GC ---
    print("[3/6] Computing GC content ...")
    windows = compute_gc(windows, args.reference)
    print(f"  GC range: [{windows['GC'].min():.4f}, {windows['GC'].max():.4f}]")
    out_cols = ["#CHR", "START", "END", "region_id", "GC"]

    # --- Step 4: Compute mappability (optional) ---
    if args.mappability_bed:
        print("[4/6] Computing mappability ...")
        windows["MAP"] = compute_mappability(
            windows,
            args.mappability_bed,
            args.genome_size,
        )
        out_cols.append("MAP")
    else:
        print("[4/6] Skipping mappability (no --mappability_bed)")

    # --- Step 5: Compute replication timing (optional) ---
    if args.repliseq:
        if args.reference_version not in REPLISEQ_REFVERS:
            print(
                f"[5/6] Skipping repli-seq (unsupported for {args.reference_version})"
            )
        else:
            print("[5/6] Computing replication timing ...")
            windows["REPLI"] = compute_repliseq(windows, standard_chroms, args)
            out_cols.append("REPLI")
    else:
        print("[5/6] Skipping repli-seq (--repliseq not set)")

    # --- Summary ---
    print_summary(windows)

    # --- Step 6: Sort and write output ---
    windows = _sort_windows(windows)
    os.makedirs(os.path.dirname(os.path.abspath(args.out_file)), exist_ok=True)
    windows[out_cols].to_csv(
        args.out_file,
        sep="\t",
        header=True,
        index=False,
        compression="gzip" if args.out_file.endswith(".gz") else None,
    )
    elapsed = time.time() - t0
    print(f"[6/6] Wrote {len(windows)} rows to {args.out_file}")
    print(f"  Total elapsed: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
