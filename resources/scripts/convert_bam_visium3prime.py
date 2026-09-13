"""Convert a Visium HD 3' BAM from <src_col> resolution to <tgt_col> resolution.

Last update: 2026-07-26

Inputs:
- argv: parquet_file, bam_in, bam_out, src_col, tgt_col
- parquet_file: the Space Ranger square_XXXum barcode mapping
Outputs:
- bam_out: the input BAM with CB retagged to the target resolution
"""

import sys

import pysam
import pandas as pd


[_, parquet_file, bam_in, bam_out, src_col, tgt_col] = sys.argv

assert src_col in {"square_002um", "square_008um", "square_016um"}, (
    f"unknown src_col: {src_col}"
)
assert tgt_col in {"square_002um", "square_008um", "square_016um"}, (
    f"unknown tgt_col: {tgt_col}"
)

# 1) read mapping
bm = pd.read_parquet(parquet_file, columns=[src_col, tgt_col])
# optional: keep only squares that are assigned to a cell (recommended)
# if "in_cell" in bm.columns:
#     bm = bm[bm["in_cell"] == 1]

print(f"#{src_col}={len(bm[src_col].unique())}")
print(f"#{tgt_col}={len(bm[tgt_col].unique())}")

src_to_tgt = dict(zip(bm[src_col].astype(str), bm[tgt_col].astype(str)))
print("mapping size:", len(src_to_tgt))

report_every = 1_000_000
kept = 0
processed = 0
missing_cb = 0

fin = pysam.AlignmentFile(bam_in, "rb")
fout = pysam.AlignmentFile(bam_out, "wb", template=fin)
for r in fin.fetch(until_eof=True):
    processed += 1
    if processed % report_every == 0:
        print(
            "processed alignments:",
            processed,
            "kept alignments:",
            kept,
            "alignments missing CB: ",
            missing_cb,
        )
    if not r.has_tag("CB"):
        missing_cb += 1
        continue

    cb = r.get_tag("CB")
    tgt = src_to_tgt.get(r.get_tag("CB"), None)
    if tgt is None:
        continue
    r.set_tag("CB", tgt, value_type="Z")
    fout.write(r)
    kept += 1
print("reads kept:", kept)

fin.close()
fout.close()

# 3) index
pysam.index(bam_out)
print("wrote", bam_out)
