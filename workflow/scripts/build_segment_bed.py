"""Stamp each configured segment with its chromosome arm and subtract the blacklist.

Last update: 2026-08-07

Inputs:
- segment_bed: configured segmentation BED4, 4th column is seg_id
- region_bed: chromosome arms BED4, 4th column is region_id
- blacklist_bed: optional BED3, subtracted from every segment
Outputs:
- aux_dir/segment.bed: BED5 carrying region_id and seg_id
"""

import logging
import os

from io_utils import read_BED
from range_utils import assign_range_to_range, trim_range_by_range
from utils import log_hist, setup_logging

snakemake_handle = snakemake  # noqa: F821
setup_logging(snakemake_handle.log[0])

segments_bed = snakemake_handle.input["segments"]
region_bed = snakemake_handle.input["region_bed"]
blacklist_bed = snakemake_handle.input.get("blacklist_bed", None)
out_bed = snakemake_handle.output["segment_bed"]

ID_COLS = ["region_id", "seg_id"]


# configured segmentation, chr-normalized (seg_id = 4th column, else CHR:START-END)
segments = read_BED(segments_bed, col_id="seg_id")
assert len(segments) > 0, f"segment_bed, no rows: {segments_bed}"
assert not segments["seg_id"].duplicated().any(), (
    f"segment_bed, duplicate seg_id: {segments_bed}"
)

# chromosome arms; every segment must lie inside exactly one of them
arms = read_BED(region_bed, col_id="region_id")
segments, na_idx = assign_range_to_range(segments, arms, "region_id", rule="contained")
assert len(na_idx) == 0, (
    f"segment_bed, {len(na_idx)} segment(s) not inside one region_bed arm: "
    f"{segments.iloc[na_idx][['#CHR', 'START', 'END']].head(3).to_dict('records')}"
)
logging.info(
    f"{len(segments)} segments over {segments['region_id'].nunique()} arms "
    f"({os.path.basename(segments_bed)})"
)

if blacklist_bed:
    bl = read_BED(blacklist_bed)
    n_raw = len(segments)
    segments = trim_range_by_range(segments, bl)
    logging.info(
        f"blacklist subtracted: {n_raw} -> {len(segments)} pieces "
        f"({os.path.basename(blacklist_bed)})"
    )

out = segments[["#CHR", "START", "END", *ID_COLS]]
out = out.sort_values(["#CHR", "START"]).reset_index(drop=True)

logging.info(
    f"segment bed: {len(out)} rows, {out['seg_id'].nunique()} seg_ids, "
    f"{out['region_id'].nunique()} region_ids"
)
out.to_csv(out_bed, sep="\t", header=False, index=False)

log_hist((out["END"] - out["START"]) / 1000.0, "segment length (kbp)")
