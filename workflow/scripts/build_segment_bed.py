"""Cut the chromosome arms at the SV extremities and subtract the blacklist.

Last update: 2026-08-17

Inputs:
- region_bed: chromosome arms BED4, 4th column is region_id
- extremity_tsv: optional headered TSV of SV breakpoints, `#CHR` and `POS0`
- blacklist_bed: optional BED3, subtracted from every segment
Outputs:
- aux_dir/segment.bed: BED5 carrying region_id and seg_id

Notes/References:
- An extremity at POS0=p cuts its arm into [START, p) and [p, END), so no bin and no bb
  spans a breakpoint and no base is lost.
- An extremity outside every arm, or inside the blacklist, is skipped.
"""

import logging
import os

from io_utils import read_BED, read_extremity_tsv
from range_utils import (
    assign_pos_to_range,
    overlaps_any_range,
    split_range_at_pos,
    trim_range_by_range,
)
from utils import log_hist, setup_logging

snakemake_handle = snakemake  # noqa: F821
setup_logging(snakemake_handle.log[0])

region_bed = snakemake_handle.input["region_bed"]
extremity_tsv = snakemake_handle.input.get("extremity_tsv", None)
blacklist_bed = snakemake_handle.input.get("blacklist_bed", None)
out_bed = snakemake_handle.output["segment_bed"]

ID_COLS = ["region_id", "seg_id"]


# chromosome arms; the segmentation is a partition of them
arms = read_BED(region_bed, col_id="region_id")
assert len(arms) > 0, f"region_bed, no rows: {region_bed}"
assert not arms["region_id"].duplicated().any(), (
    f"region_bed, duplicate region_id: {region_bed}"
)
logging.info(f"{len(arms)} arms ({os.path.basename(region_bed)})")

blacklist = read_BED(blacklist_bed) if blacklist_bed else None

if extremity_tsv:
    ext = read_extremity_tsv(extremity_tsv)
    n_raw = len(ext)
    ext, off_arm = assign_pos_to_range(ext, arms, "region_id", dropna=True)
    logging.info(f"extremities outside every arm, skipped: {len(off_arm)}/{n_raw}")
    if blacklist is not None:
        hit = overlaps_any_range(ext, blacklist)
        logging.info(
            f"extremities inside the blacklist, skipped: {hit.sum()}/{len(hit)}"
        )
        ext = ext[~hit]
    segments = split_range_at_pos(arms, ext)
    logging.info(
        f"{len(ext)}/{n_raw} extremities cut {len(arms)} arms into {len(segments)} "
        f"segments ({os.path.basename(extremity_tsv)})"
    )
else:
    segments = arms
    logging.info("no extremity_tsv; one segment per arm")

segments["seg_id"] = (
    segments["region_id"]
    + "#"
    + segments["START"].astype(str)
    + "-"
    + segments["END"].astype(str)
)

if blacklist is not None:
    n_raw = len(segments)
    segments = trim_range_by_range(segments, blacklist)
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
