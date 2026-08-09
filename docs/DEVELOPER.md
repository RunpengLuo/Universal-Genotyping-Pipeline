# Developer guide: vocabulary and naming

The words this codebase uses for genomic units, matrix axes and function names. One word
per concept. Read this before adding a module; a new synonym is a bug.

User-facing names (config keys, output columns, output filenames, rule names) are frozen
and may keep older wording - see [Legacy surface names](#legacy-surface-names).

## Table of Contents

- [Genomic units](#genomic-units)
- [Matrix axes](#matrix-axes)
- [Supporting words](#supporting-words)
- [Identifiers](#identifiers)
- [Retired synonyms](#retired-synonyms)
- [Legacy surface names](#legacy-surface-names)
- [Function naming](#function-naming)

## Genomic units

Four nested levels, plus two entities that are assigned into them.

```
region  (chromosome arm, region_id)      <- region_bed column 4
  segment (seg_id from segment_bed)      <- build_segment_bed; the hard bb boundary
    bin     (fixed tile, bin_id)         <- window_size tiles; mosdepth + rd_correct
      bb      (merged bin, bb_id)        <- build_adaptive_bins; the output unit
SNP     (allele axis entry, POS0)        <- assigned into a bin and a bb
gene / exon  (GTF entities, feature_id)
```

| Unit | Id | Defined by | Meaning |
|---|---|---|---|
| region | `region_id` | `region_bed` column 4 | Chromosome arm. Carried for RDR/QC; a bb belongs to exactly one. |
| segment | `seg_id` | `build_segment_bed` | The 4th column of `segment_bed`. Binning clusters by it and never merges across it. With `segment_bed == region_bed` there is one segment per arm. |
| bin | `bin_id` | `build_window_bed` | Fixed `window_size` tile (1 kb). What mosdepth counts and `rd_correct` bias-corrects. Internal; never written to an output. |
| bb | `bb_id` | `build_adaptive_bins` | Consecutive bins merged until every tumor observation meets `min_snp_reads`. The output unit: features of `bb.tsv.gz` and every `bb.*.npz`. |
| SNP | `POS0` | genotyping + phasing | Phased het SNP, the allele axis. Assigned into its containing bin and bb. |
| gene / exon | `feature_id` | the GTF | Annotation only. `feature_type` is exon > intron > intergenic. |

A **bb is a bin too, just larger**. Say "fixed bin" and "merged bin" only where the
distinction matters; otherwise `bin` and `bb` carry it.

## Matrix axes

Every matrix in the pipeline is `(n_features, n_observations)`. Never say "row" or
"column" for them.

| Word | Meaning | Instances |
|---|---|---|
| `feature` | the genomic axis entry | bin, bb, SNP, gene, exon, peak |
| `observation` | the measurement axis entry | cell, spot, dataset/replicate, bulk sample |

`feature_id` / `feature_type` are the GTF-derived output columns: a gene is one kind of
feature. Matplotlib subplot geometry (`nrows`, `ncols`, `axes[ri, si]`) is figure layout,
not a data axis, and keeps its own words.

## Supporting words

| Word | Meaning |
|---|---|
| `pos` | a 0-based coordinate (`POS0`) |
| `range` | any `[START, END)` pair, 0-based half-open, on either the query or the reference side |
| `cluster` | a set of items processed together: summation groups, phase clusters (VCF `PS`), non-overlapping range clusters |
| `span` | bp length of a unit (a measure, not a unit) |

## Coordinate invariant

**Every pos is 0-based. Every range is 0-based, left-closed, right-open `[START, END)`.**
`range_utils` asserts both on every call: `pos_col="POS"` is rejected outright, and a
reference or query range with `END <= START` fails with a named message.

A 1-based value exists only at these boundaries, and is converted on the spot:

| Site | Why |
|---|---|
| `read_VCF` | VCF `POS` is 1-based; it derives `POS0 = POS - 1` |
| `read_bcftools_pileup_counts` | `bcftools query` emits 1-based `POS` |
| `read_gtf` | GTF is 1-based inclusive; converted to `[START, END)` on ingest |
| `parse_genetic_map` | published genetic maps are 1-based |
| `interp_cM_between_bbs` | interpolates against the 1-based map, so it converts the bb hull with `START + 1` / `END` |
| the `POS` column of `snps.tsv.gz` | output kept 1-based for VCF compatibility |

Everywhere else, `POS` is data to carry, never a coordinate to compute with.

## Identifiers

Three suffixes, three meanings:

| Suffix | Meaning | Examples |
|---|---|---|
| `_id` | a unit's identifier | `region_id`, `seg_id`, `bin_id`, `bb_id`, `snp_id`, `feature_id`, `cluster_ids` |
| `_df_idx` | a row index into a DataFrame, and into any matrix aligned to it | `_orig_df_idx`, `RAW_SNP_DF_IDX` |
| `_idx` | a positional index local to one function, into neither | `hit_idx`, `ideal_idx`, `gc_bin_idx` |

An id survives a reindex; a `_df_idx` does not. Anything crossing a function or column
boundary must say which of the two it is.

## Retired synonyms

Do not reintroduce these:

| Retired | Use instead |
|---|---|
| window | bin (the fixed tile) |
| block | bb, cluster, or span, depending on which of its five old senses is meant |
| interval | range |
| grid | the unit's own name, or drop the word |
| group | cluster |
| unit | the unit's own name |
| tier | cluster |
| row / column / col | feature / observation |
| chunk (genomic) | segment. The pandas `chunksize` read batch keeps the word. |
| rep / replicate (in an identifier) | dataset. `dataset_id` is the sample-sheet key and the `{dataset_id}` wildcard; prose may still say "replicate". |
| arm | region (prose may still say "chromosome arm") |
| segmentation (meaning binning) | binning; `seg_id` is the unrelated segment level |

## Legacy surface names

These predate the vocabulary and are frozen because they are user-facing. The code around
them uses the canonical word.

| Surface name | Canonical concept |
|---|---|
| `window_size`, `window_bed`, `params_build_windows` (config keys) | fixed bin |
| `aux/windows.bed.gz`, `windows.3col.bed.gz`, `window.dp.npz`, `window.tsv.gz` (outputs) | fixed bin |
| `build_windows.smk`, `build_window_bed`, `window_bed_to_3bed` (rules) | fixed bin |
| `BLOCKSIZE` (output column), `max_blocksize` (config key) | bb span |
| `min_snp_per_bin`, `min_snp_reads` (config keys) | thresholds on a **bb** |
| `bb.tsv.gz`, `bb_dir`, `bb_id` | the merged bin (HATCHet's `.bb` convention) |
| `PS` (VCF tag and column) | phase cluster |
| `REP_ID`, `RDR_BASE_REP_ID` (output columns) | `dataset_id`. See the note below. |

> [!NOTE]
> The window BED **is** the fixed-bin BED. `rd_correct` reads `input["window_bed"]` into
> a local named `bin_df`.

> [!IMPORTANT]
> **Future goal: rename `REP_ID` -> `DATASET_ID` in the output files.**
> Every internal identifier now says `dataset`; `REP_ID` and `RDR_BASE_REP_ID` survive
> only as on-disk column names, in `sample_ids.tsv`, `barcodes.full.tsv.gz`, and the
> `{BARCODE}_{REP_ID}` values of `barcodes.tsv.gz` (see `docs/reference.md`). Renaming
> them is a breaking output change for HATCHet / Copy-typing / CalicoST, so it needs its
> own release. Until then a `"REP_ID"` string literal next to a `dataset_id` variable is
> expected, not an oversight.

## Function naming

`<verb>_<object>[_to_<target>]`, verb first, target unit named explicitly.

| Verb | Contract |
|---|---|
| `read_*` | parse a file (`io_utils` only) |
| `build_*` | construct a unit set, id map, or descriptor; no matrix math |
| `assign_*` | attach a target-unit id to each query feature; no sums. Returns `(annotated_qry, na_idx)`, never mutates the query, and takes `fillna=` / `dropna=` (mutually exclusive). `na_idx` holds the positional indices of the unassigned rows, so a parallel matrix subsets the same way |
| `map_*` | reindex a matrix onto a target feature axis, 0-filling the absent features; no sums, the axis length changes but no value does |
| `sum_*` | collapse a matrix along one axis using a cluster-id array |
| `count_*` | tally events into units |
| `aggregate_*` | collapse with weighting or a non-trivial reduction |
| `interp_*` | derive a value or bound between known points (`interp_pos_ranges`, `interp_cM_between_bbs`) |
| `compute_*` | derive a numeric quantity from already-aggregated data |
| `estimate_*` | statistical estimate |
| `mask_*` | return a boolean selection |

Two kernels sit under everything else; the rest are callers or wrappers.

| Kernel | Module | Wrappers |
|---|---|---|
| `cluster_sum` | `matrix_utils` | `sum_features_to_bbs` (feature axis), `sum_observations_to_pseudobulk` (observation axis) |
| `_merge_bins_to_bbs` (numba) | `aggregation_utils` | `build_adaptive_bins` |
| `_searchsorted_assign` | `range_utils` | `assign_pos_to_range`, `assign_pos_to_range_ovlp`, `overlaps_any_range`, and `assign_range_to_range` under `rule="midpoint"` / `"contained"` |

`assign_range_to_range` takes `rule=`: `max_overlap` (largest shared span, the only rule
that walks queries one at a time), `midpoint`, or `contained` (an id only when `START` and
`END - 1` land in the same reference **id** - two rows sharing one `seg_id` count as one,
so a window over a blacklist hole is contained, not straddling).

Worked examples:

```python
assign_pos_to_range(snps, bins, ref_id="bin_id", dropna=True)  # -> (snps.bin_id, na_idx)
assign_range_to_range(win, seg, "seg_id", rule="contained")    # -> (win.seg_id, na_idx)
build_adaptive_bins(bins, snps, ...)          # -> bbs, snps.bb_id
aggregate_bin_depth_to_bbs(...)               # length-weighted depth
sum_features_to_bbs(X, bb_ids, n_bbs)         # cluster_sum(axis=0)
sum_observations_to_pseudobulk(X, ids, n)     # cluster_sum(axis=1)
sum_atac_fragments_to_bins(...)              # tally, not a matrix collapse
```
