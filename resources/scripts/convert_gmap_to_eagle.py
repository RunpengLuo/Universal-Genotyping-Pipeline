#!/usr/bin/env python3
"""Convert a SHAPEIT per-chromosome genetic map into Eagle's single-table format.

Last update: 2026-07-26

Inputs:
- argv: see parse_args for the input maps and output path
Outputs:
- gmap: the Eagle table, chr POS rate cM
"""

import argparse
import os
import sys

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert per-chromosome SHAPEIT genetic maps to Eagle2 format."
    )
    parser.add_argument(
        "input_dir", help="Directory with per-chr gmap files (chr{N}.*.gmap.gz)"
    )
    parser.add_argument("output_file", help="Output gzipped Eagle-format genetic map")
    parser.add_argument(
        "--pattern",
        default="{chrom}.t2t.scaled.gmap.gz",
        help="Filename pattern with {chrom} placeholder (default: {chrom}.t2t.scaled.gmap.gz)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    chroms = [f"chr{c}" for c in list(range(1, 23)) + ["X"]]
    frames = []

    for chrom in chroms:
        fname = args.pattern.format(chrom=chrom)
        fpath = os.path.join(args.input_dir, fname)
        if not os.path.isfile(fpath):
            print(f"WARNING: skipping {fpath} (not found)", file=sys.stderr)
            continue

        df = pd.read_csv(fpath, sep="\t")
        chr_num = chrom.replace("chr", "")
        df.insert(0, "chr", chr_num)
        df.columns = ["chr", "position", "COMBINED_rate(cM/Mb)", "Genetic_Map(cM)"]
        frames.append(df)
        print(f"  {chrom}: {len(df)} entries")

    merged = pd.concat(frames, ignore_index=True)
    merged.to_csv(args.output_file, sep=" ", index=False, compression="gzip")
    print(f"wrote {len(merged)} entries to {args.output_file}")


if __name__ == "__main__":
    main()
