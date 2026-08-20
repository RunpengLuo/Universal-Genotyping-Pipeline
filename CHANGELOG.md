# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0b2] - 2026-08-13

One explicit segmentation and one shared bin grid in every mode, references without a
`chr` prefix, and a flat `allele_dir`. Output paths, column names and bb boundaries move.

### Added

#### Config
- `extremity_tsv`: a TSV of SV breakpoints (`#CHR`, `POS0`) that no bin or bb may cross.
- `extremity_tsv` unset leaves the arms uncut, one segment per chromosome arm.
- `build_segment_bed` cuts the `region_bed` arms at each breakpoint and subtracts the
  blacklist; a breakpoint off-arm or in the blacklist is skipped.
- `extremity_tsv` ignores a pre-built `window_bed` and re-tiles from the cut segments.
- `resources/templates/extremity.tsv`: an example breakpoint file.
- `dataset_ids` (default `[]`, all): restrict a run to a subset of the `sample_id`'s
  datasets. A multiome pair shares one `dataset_id`, so naming it keeps both records;
  a named `dataset_id` with no record on the run's `reference_version` is an error.
- Blacklist pieces of a segment keep its `seg_id`, so holes never bound a bin.
- `species` (required, `human` | `mouse`) sets the sex-chromosome numbering.
- `REFVERS_ALIAS` folds spellings: `GRCh38` -> `hg38`, `T2T-CHM13v2.0` -> `chm13v2`.
- An alias now also enables build-gated steps such as Repli-seq correction.

#### Workflow
- References whose contigs lack a `chr` prefix, the style read from `genome_size`.
- The prefix is converted only at tool boundaries; internal frames stay chr-prefixed.
- Repli-seq correction on chm13v2, lifted from hg19 by the UCSC `hg19ToHs1` chain.
- Bulk multi-SNP diagnostics at `bb_dir/multi_snp/bulk/`, with depth and RDR.
- `rd_correct` also writes `pileup_dir/bulk/window.raw.dp.npz`, the uncorrected,
  unmasked mosdepth depth on the same axes as `window.dp.npz`.

#### Development
- `docs/DEVELOPER.md`: the vocabulary, coordinate conventions and module map.
- Unit tests for the range primitives, refver folding, contig naming, IO readers and the
  sample-sheet loader.
- `ruff` and `snakefmt --check workflow/` run in CI.

### Changed

#### Output layout (breaking)
- `allele_dir` is flat and identical in every mode; its subdirs are gone.
- Single-cell allele matrices are one union over every assay, not one per assay.
- The single-cell column key is `{barcode}_{dataset_id}_{assay_type}`; h5ad `obs_names`
  follow.
- `bb_dir` keeps its per-assay subdirs, sliced out of that union.
- `sample_ids.tsv` columns are the record keys: `sample_id`, `dataset_id`,
  `rdr_base_dataset_id`.
- The h5ad `obs` follows; only the derived `SAMPLE` keeps an uppercase name.
- Observation order is fixed at parse time: assay, normal before tumor, `dataset_id`.
- Bulk corrected depth is one `pileup_dir/bulk/window.dp.npz` over every bulk dataset.
- `pileup_dir/{assay}/window.tsv.gz` is gone; the row axis is the window BED.
- One `depth_statistics.tsv` and one `qc/rd_correction.bulk.pdf` per run.
- Single-cell multi-SNP diagnostics move to `bb_dir/multi_snp/{assay_type}/`, out of the
  sweep.
- They are written by `write_bb_file`, so the bb schema replaces the `multi_id` frame.
- `copytyping_preprocess` writes `bb.tsv.gz` (bb schema plus `#feature`), not
  `cnv_segments.tsv`.
- `barcodes.full.tsv.gz` is dropped; the dataset splits off `barcodes.tsv.gz`.
- `allele_dir/{assay}/unique_snp_ids.npy` is dropped, no consumer.

#### Binning and read-depth correction (values move)
- Single-cell binning uses the window BED as fixed bins, the grid bulk bins on.
- `build_window_bed` therefore runs in every mode.
- Single-cell bb boundaries follow 1 kb windows rather than SNP density.
- A segment holding no het SNP now yields bbs with `#SNPS = 0`.
- scATAC `Xcount` drops: fragments count through windows, never in a blacklist hole.
- `seg_id` is `{region_id}#{START}-{END}`, no longer `{region_id}#{k}`.
- Single-cell binning clusters on `seg_id`, so its bbs stop at segment bounds.
- Non-bulk SNP ranges are bounded at segment edges, no longer spanning blacklist holes.
- Zero raw depth corrects to `0.0`; NaN now means only "correction undefined".
- The mappability floor moved into `rd_correct`, so it applies without GC correction too.
- `correct_readcount_lowess` no longer extrapolates outside the fitted covariate range.
- **Breaking**: `params_count_reads.gc_correct_method` is renamed `rd_correct_method`.
- `rd_correct_method: median` values are unchanged; `lowess` values shift.
- `rd_correct` runs once per bulk run instead of once per assay.
- Correction NaN is masked per column, not dropped for the whole assay.
- `min_snp_reads` thresholds tumor columns only in single-cell, as bulk already did.
- A shallow single-cell normal no longer coarsens binning (synthetic case: 16 bbs, not 2).
- A run with no tumor warns and falls back to `min_snp_per_bin` alone.
- Single-cell multi-SNP groups are built from window-assigned SNPs, as bulk already did.
- Per-SNP `START`/`END` are computed once over the shared grid, so group bounds move.
- Default `min_snp_reads` sweeps eight values `100..10000`; `min_snp_per_bin` is `1`.
- Only bulk gets the GC/MAP/REPLI covariates in its window BED.
- The Repli-seq fetch now also requires `params_count_reads.rt_correct`.

#### Sample sheet (breaking)
- `reference_version` is required on every record.
- A run keeps only records whose build folds to the config's `reference_version`.
- An unset config `reference_version` is now an error, not a warning.
- The TSV sheet is a flat JSON schema: one `files.<key>` column per input.
- TSV sheets can now carry remote URLs and name single-cell files explicitly.
- `chromosomes` is validated against `genome_size`; a missing chromosome is an error.
- Every sample-file input is checked at DAG build; a URL is never fetched.
- `sample_id` and `dataset_id` must match `[A-Za-z0-9_-]+`, checked as the sheet loads.
- It replaces the `dataset_id`-only check in `parse_records`, which took non-ASCII
  letters and never saw `sample_id`.

#### Genotyping
- `apply_clonal_loh_hmm` (default `false`, bulk only, in `params_genotype_snps`)
  re-genotypes a high-purity tumor's calls with a latent clonal-LOH chain, rescuing
  germline hets that LOH collapsed onto one allele. Every genotyped dataset must be
  `tumor`. Off, the calls symlink through unchanged.
- On, the caller keeps every callable panel site (`--keep-alts`, no `--variants-only`, no
  `GT="alt"`), so a `0/0` site survives and stays phaseable once re-called.
- The chain is `script_utils/genotype_loh_hmm.py`; one `region_bed` arm is one chain.
  Parameters and guidance: `docs/reference.md`.
- It writes `aux_dir/clonal_loh_hmm.{segments,params}.tsv`, header-only when it did not run.
- `params_annotate_snps` becomes `params_genotype_snps`.
- `min_hom_dp` is gone. `min_dp` is one depth floor for both modes and now gates het calls
  too, so a single-cell site below it is a no-call.
- `annotate_snps_pseudobulk.py` becomes `post_genotype_snps.py`, one branch per mode.
- `post_genotype_snps` carries `ID`, `QUAL` and `FILTER` through from its input VCF.
- The genotyping QC plot colours het blue, hom red, no-call grey.
- The bulk caller writes `snp_dir/raw/chr{chrname}.vcf.gz`; `snp_dir/chr{chrname}.vcf.gz`
  is post-processing's output, so phasing reads the same path as before.
- `bcftools call` is always constrained to the panel's REF/ALT (`--constrain alleles`), so
  a somatic allele can never become the called ALT. The panel REF must match `reference`.
- `genotype_snps_bulk` pipes `call` into the filtering `view` instead of staging an
  unfiltered VCF, under `set -euo pipefail`.

#### Genotyping (breaking)
- Bulk genotyping reads `snp_panel` directly, via `bcftools mpileup -T`, which takes
  CHROM/POS only: REF comes from `reference`, ALT from the reads.
- `snp_panel` is required in both genotyping modes and must be a bgzipped, indexed VCF.
  `.bcf` is refused at parse time (samtools/bcftools#690).
- `genotype_snps_bulk` and `pileup_snps_bulk_bcftools` pass `--regions` in every
  `remote_mode`, not only `stream`.
- Every bcftools option in `workflow/rules/` is spelled long.

#### Performance (same outputs)
- Bulk pileup splits per chromosome and is joined by `merge_pileup_counts`; the
  concatenation is byte-identical to the single-job counts file.
- The clonal-LOH chain's recursions are numba kernels.
- `--threads` moves off the single-threaded `bcftools mpileup` onto the compressing steps;
  `threads.genotype` 4 -> 2 and `threads.pileup` 8 -> 2.
- Genotyping counts its records with `bcftools query --format '\n'`.
- `bcftools call --variants-only` when the HMM is off, which only shrinks the temporary
  unfiltered VCF.

#### Internals
- One word per concept: `region` > `segment` > `bin` > `bb`, `feature` x `observation`.
- `rep` becomes `dataset_id` in every identifier.
- Five range-assignment loops collapse into `range_utils.py` under one contract.
- 0-based half-open coordinates are enforced rather than assumed.
- Rules read `parse_workflow` globals; only `params_*` and `threads` stay `config[...]`.
- No module uses `import *`, so an undefined name now fails CI.
- Every mask in `combine_counts` logs as `label: n/total (%unit=0.xxx) (xxx.xxx Mbps)`.
- The GTF is parsed once per run.
- SNPs are assigned to fixed bins once, not once per sweep point.

#### Resources
- Bundled `region_bed` `region_id` is an arm label (`chr1p`), no longer `CHR:START-END`.
- Pre-built window BEDs carry those labels and `seg_id` `{region_id}#{START}-{END}`.

### Fixed
- **Output change**: a bedGraph midpoint in a window gap was credited to its predecessor.
- **Output change**: duplicate bedGraph rows in one window overwrote rather than summed
  `REPLI`.
- **Output change**: per-SNP ranges compared 1-based `POS` to 0-based bounds, off by one.
- **Output change**: non-bulk `snps.tsv.gz` now carries `PS`; binning was one phase
  cluster.
- A fixed bin with a null cluster key was folded into bb 0.
- Per-bin `PS` labels are filled both ways; a null cluster key is an assertion.
- A SNP-free bb produced NaN switch probabilities and propagated them to the next bb.
- On a bare-contig genome `build_window_bed` dropped every bin as off-segment.
- `mappability_bed` ranges reached bedtools unrenamed.
- An Ensembl GTF annotated every SNP `intergenic`, silently degrading
  `gene_aware_binning`.
- Contigs are no longer renamed on disk.
- `chromosomes: [chr22]` produced `phase/chrchr22.vcf.gz`.
- `parse_genetic_map` no longer guesses which numeric contig is X.
- A chromosome still missing after relabeling is an error naming it.
- The bin BED derives its contigs from `genome_size`, not a hardcoded `chr1..chr22`.
- A JSON `files` entry set to `null` became the literal path `"None"`.
- A record missing a required key raised a bare `KeyError`.
- UCSC downloads use `hgdownload.soe.ucsc.edu`; the `cse` host failed TLS verification.
- `phase_and_concat_nonbulk` used `scanpy` without importing it.

### Removed
- **Breaking**: SV-breakpoint BEDPE support; `extremity_tsv` lists the breakpoints
  instead.
- **Breaking**: `rdr_outlier_quantile` and its clipping, which discarded focal
  amplifications.
- **Breaking**: the genotyping QC plot, `qc_genotype_snps` and `quality_control.smk`.
- **Breaking**: `pseudobulk_snp_statistics.tsv` and its single-cell `report()` entry.
- **Breaking**: `snp_targets`; `snp_panel` is the only SNP-site input in every mode.
- `resources/scripts/build_snp_targets.sh`, and the `target_positions/` output of the
  1kGP and mouse MGP panel builders.
- `verify_window_bed` and the `window_bed.checked` gate; the build cannot cross a segment.
- `resources/scripts/validate_sample_file.py`; the DAG build runs the same validation.
- The `tests/data/` end-to-end cases (COLO829, HCC1395) and their CI dry-run step.

## [0.1.0b1] - 2026-07-26

Initial pre-release.

### Added

#### Workflow
- Snakemake pipeline for SNP genotyping, phasing, and allele counting across three
  modes (`bulk_genotyping`, `single_cell_genotyping`, `copytyping_preprocess`), with
  per-rule conda environments, logs, and benchmarks.
- JSON sample file (legacy TSV supported); remote `http(s)` inputs downloaded
  (`remote_mode: storage`) or streamed per chromosome (`stream`, bulk only).
- Per-mode tutorials, config/output reference, sample-sheet spec, and resource catalog.

#### Common rules
- Phasing via Eagle2, SHAPEIT5, or LongPhase; a supplied `het_snp_vcf` short-circuits
  genotyping (and phasing when already phased).
- Multi-replicate allele and depth consolidation onto one shared SNP and bb set.

#### `bulk_genotyping`
- Bulk assays (WGS/WGS-lr/WES) on one shared bin/bb set, WES handled like WGS.
- bcftools genotyping from matched normals and het-SNP pileup; mosdepth depth with
  GC/mappability/replication-timing correction; adaptive binning to RDR/BAF matrices
  for HATCHet3.

#### `single_cell_genotyping`
- scRNA, scATAC, VISIUM/VISIUM3prime assays, including 10x Epi Multiome pairs.
- cellsnp-lite pseudobulk genotyping and per-cell pileup to allele and native-count
  matrices for CalicoST.

#### `copytyping_preprocess`
- Aggregate single-cell / spatial allele and native counts onto a pre-computed phased
  VCF and bb set (no genotyping or phasing) for CalicoST.
