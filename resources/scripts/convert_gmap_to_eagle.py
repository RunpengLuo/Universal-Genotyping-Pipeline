#!/usr/bin/env python3
"""Convert per-chromosome SHAPEIT genetic maps to Eagle2 format.

Eagle2 expects a single gzipped file with columns:
  chr  position  COMBINED_rate(cM/Mb)  Genetic_Map(cM)

SHAPEIT per-chromosome files (e.g., chr1.t2t.scaled.gmap.gz) have columns:
  pos  cM/Mb  cM

Usage:
  python convert_gmap_to_eagle.py <input_dir> <output_file>

Example:
  python convert_gmap_to_eagle.py \\
    /path/to/t2t_native_scaled_maps \\
    /path/to/genetic_map_chm13v2_withX.txt.gz
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
