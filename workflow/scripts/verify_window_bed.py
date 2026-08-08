"""Check a supplied window BED against the segment BED.

Runpeng Luo
Last update: 2026-08-07

Runs only when config `window_bed` is set, i.e. the windows are read rather than tiled
from the segments. A window that crosses a segment bound would let binning merge a bb
across a segment, so it is a hard error; the run must either fix the file or unset
`window_bed` to rebuild. The window BED's own region_id/seg_id are assumed to be the
segment_bed ids and are not rewritten.

Dependencies:
  pandas; io_utils.read_BED, range_utils.assign_range_to_range,
  range_utils.overlaps_any_range, utils.setup_logging, utils.add_chr_prefix.

Inputs
  window_bed: window BED (headered, #CHR START END region_id seg_id GC [MAP] [REPLI]).
  segment_bed: BED5 from build_segment_bed.
Outputs:
  checked: one-line summary; its presence gates the read-depth rules.
"""

import logging

import pandas as pd

from io_utils import read_BED
from range_utils import assign_range_to_range, overlaps_any_range
from utils import add_chr_prefix, setup_logging

snakemake_handle = snakemake  # noqa: F821
setup_logging(snakemake_handle.log[0])

window_bed = snakemake_handle.input["window_bed"]
segment_bed = snakemake_handle.input["segment_bed"]
chroms = list(snakemake_handle.params["chroms"])
out_checked = snakemake_handle.output["checked"]

windows = pd.read_table(window_bed, sep="\t", dtype={"#CHR": str})
assert "#CHR" in windows.columns, (
    f"window_bed, missing the `#CHR ...` header row: {window_bed}"
)
windows["#CHR"] = add_chr_prefix(windows["#CHR"])
n_raw = len(windows)
windows = windows[windows["#CHR"].isin(chroms)].reset_index(drop=True)
logging.info(f"{len(windows)}/{n_raw} windows on {len(chroms)} chromosomes")
assert len(windows) > 0, f"window_bed, no window on {chroms[:3]}...: {window_bed}"

# a window is contained when both of its ends land in the same seg_id; a window lying in
# no segment at all (a blacklist gap) is not a crossing
segments = read_BED(segment_bed)[["#CHR", "START", "END", "seg_id"]]
windows, _ = assign_range_to_range(
    windows, segments, "seg_id", rule="contained", out_col="_seg"
)
windows["_end0"] = windows["END"] - 1
outside = ~overlaps_any_range(windows, segments, pos_col="START") & ~overlaps_any_range(
    windows, segments, pos_col="_end0"
)
crossing = windows["_seg"].isna().to_numpy() & ~outside

assert not crossing.any(), (
    f"window_bed, {int(crossing.sum())} window(s) cross a segment bound: "
    f"{windows.loc[crossing, ['#CHR', 'START', 'END']].head(3).to_dict('records')}"
)
if outside.any():
    logging.warning(
        f"WARN: {int(outside.sum())} window(s) lie in no segment (blacklist gaps); "
        "they are counted and binned under their own ids"
    )

summary = (
    f"windows={len(windows)}\tsegments={len(segments)}\t"
    f"crossing=0\toutside={int(outside.sum())}\n"
)
logging.info(summary.replace("\t", " "))
with open(out_checked, "w") as fh:
    fh.write(summary)
