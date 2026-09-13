#!/usr/bin/env python3
"""Build an mm10 genetic map from the published mouse recombination rates.

Last update: 2026-08-06

Inputs:
- argv: see parse_args for the source table and output path
Outputs:
- gmap: #CHR POS cM, one row per map position
References:
- resources/README.md: how the mouse panel and map are built
"""

import argparse
import gzip
import os

import numpy as np
import pandas as pd

MOUSE_CHROMS = [str(c) for c in range(1, 20)] + ["X"]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build mm10 Eagle2 + SHAPEIT5 genetic maps from CoxMapV3."
    )
    parser.add_argument("input_csv", help="CoxMaps_rev_build38.csv (CoxMapV3)")
    parser.add_argument("out_dir", help="output directory (eagle2/ + shapeit5/)")
    return parser.parse_args()


def main():
    args = parse_args()
    eagle_dir = os.path.join(args.out_dir, "eagle2")
    shapeit_dir = os.path.join(args.out_dir, "shapeit5")
    os.makedirs(eagle_dir, exist_ok=True)
    os.makedirs(shapeit_dir, exist_ok=True)

    df = pd.read_csv(
        args.input_csv,
        usecols=["chr_b38", "bp_b38", "fem_cM", "ave_cM"],
        dtype={"chr_b38": str},
    )
    df = df[df["bp_b38"] > 0]

    eagle_frames = []
    for chrom in MOUSE_CHROMS:
        sub = df[df["chr_b38"] == chrom].sort_values("bp_b38")
        if sub.empty:
            print(f"WARNING: no markers for chr{chrom}")
            continue
        pos = sub["bp_b38"].to_numpy(dtype=np.int64)
        cm_col = "fem_cM" if chrom == "X" else "ave_cM"
        cm = sub[cm_col].to_numpy(dtype=np.float64)
        # drop markers whose cM dips below the running max (CoxMapV3 chr10/chr14 inversions)
        keep = np.ones(len(cm), dtype=bool)
        keep[1:] = cm[1:] >= np.maximum.accumulate(cm)[:-1]
        pos, cm = pos[keep], cm[keep]
        rate = np.zeros_like(cm)
        dpos = np.diff(pos)
        rate[:-1] = np.where(dpos > 0, np.diff(cm) / dpos * 1e6, 0.0)

        shapeit = pd.DataFrame({"pos": pos, "chr": f"chr{chrom}", "cM": cm})
        with gzip.open(
            os.path.join(shapeit_dir, f"chr{chrom}.mm10.gmap.gz"), "wt"
        ) as fh:
            shapeit.to_csv(fh, sep="\t", index=False, float_format="%.6f")

        eagle_frames.append(
            pd.DataFrame(
                {
                    "chr": chrom,
                    "position": pos,
                    "COMBINED_rate(cM/Mb)": rate,
                    "Genetic_Map(cM)": cm,
                }
            )
        )
        print(f"  chr{chrom}: {len(pos)} markers ({int((~keep).sum())} dropped)")

    eagle = pd.concat(eagle_frames, ignore_index=True)
    out_eagle = os.path.join(eagle_dir, "genetic_map_mm10_withX.txt.gz")
    with gzip.open(out_eagle, "wt") as fh:
        eagle.to_csv(fh, sep=" ", index=False, float_format="%.6f")
    print(f"wrote {len(eagle)} markers to {out_eagle}")


if __name__ == "__main__":
    main()
