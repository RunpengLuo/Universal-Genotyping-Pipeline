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
  segment (breakpoint chunk, seg_id)     <- build_segment_bed; the hard bb boundary
    bin     (fixed tile, bin_id)         <- window_size tiles; mosdepth + rd_correct
      bb      (merged bin, bb_id)        <- build_adaptive_bins; the output unit
SNP     (allele axis entry, POS0)        <- assigned into a bin and a bb
gene / exon  (GTF entities, feature_id)
```

| Unit | Id | Defined by | Meaning |
|---|---|---|---|
| region | `region_id` | `region_bed` column 4 | Chromosome arm. Carried for RDR/QC; a bb belongs to exactly one. |
| segment | `seg_id` | `build_segment_bed` | `{region_id}#{k}`, the arm split at `breakpoint_bedpe` cuts. Binning clusters by it and never merges across it. With no BEDPE there is one segment per arm. |
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
| `range` | any `[START, END)` pair, 0-based half-open, on either the query or the reference side |
| `cluster` | a set of items processed together: summation groups, phase clusters (VCF `PS`), non-overlapping range clusters |
| `span` | bp length of a unit (a measure, not a unit) |

## Identifiers

Three suffixes, three meanings:

| Suffix | Meaning | Examples |
|---|---|---|
| `_id` | a unit's identifier | `region_id`, `seg_id`, `bin_id`, `bb_id`, `snp_id`, `feature_id`, `cluster_ids` |
| `_df_idx` | a row index into a DataFrame, and into any matrix aligned to it | `_orig_df_idx`, `feature_df_idx`, `RAW_SNP_DF_IDX` |
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
| arm | region (prose may still say "chromosome arm") |
| segmentation (meaning binning) | binning; `seg_id` is the unrelated breakpoint level |

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
| `feature_df_idx` (temp column in `assign_features_to_ranges`) | a positional index, not an id; `feature_id` is taken by the GTF gene |

> [!NOTE]
> The window BED **is** the fixed-bin BED. `rd_correct` reads `input["window_bed"]` into
> a local named `bin_df`.

## Function naming

`<verb>_<object>[_to_<target>]`, verb first, target unit named explicitly.

| Verb | Contract |
|---|---|
| `read_*` | parse a file (`io_utils` only) |
| `build_*` | construct a unit set, id map, or descriptor; no matrix math |
| `assign_*` | attach a target-unit id to each query feature; no sums |
| `sum_*` | collapse a matrix along one axis using a cluster-id array |
| `count_*` | tally events into units |
| `aggregate_*` | collapse with weighting or a non-trivial reduction |
| `compute_*` | derive a numeric quantity from already-aggregated data |
| `estimate_*` | statistical estimate |
| `mask_*` | return a boolean selection |

Two kernels sit under everything else; the rest are callers or wrappers.

| Kernel | Module | Wrappers |
|---|---|---|
| `cluster_sum` | `matrix_utils` | `sum_features_to_bbs` (feature axis), `sum_observations_to_pseudobulk` (observation axis) |
| `_merge_bins_to_bbs` (numba) | `aggregation_utils` | `build_adaptive_bins` |

Worked examples:

```python
assign_snps_to_bins(snps, bins, tot_mtx)      # -> snps.bin_id
build_adaptive_bins(bins, snps, ...)          # -> bbs, snps.bb_id
aggregate_bin_depth_to_bbs(...)               # length-weighted depth
sum_features_to_bbs(X, bb_ids, n_bbs)         # cluster_sum(axis=0)
sum_observations_to_pseudobulk(X, ids, n)     # cluster_sum(axis=1)
count_atac_fragments_to_bbs(...)              # tally, not a matrix collapse
```
