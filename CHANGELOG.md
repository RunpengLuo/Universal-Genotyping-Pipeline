# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `REFVERS_ALIAS` in `config/const.py`: accepted spellings of each reference version
  (`GRCh38` -> `hg38`, `T2T-CHM13v2.0` -> `chm13v2`, ...), matched exactly and
  case-insensitively. Both the config value and every record's `reference_version` fold
  through it.
- `validate_sample_file.py --reference-version`; without it, every genome build in the file
  is validated on its own.
- `species` config key (required; `human` | `mouse`, shipped as `human`). An unset value is
  an error, an unrecognized one warns and is used verbatim. Sex-chromosome numbering is a
  property of the species, not of the reference version.
- Support for references that name contigs without a `chr` prefix (Ensembl, b37). The style
  is read from `genome_size`; region strings, `windows.3col.bed.gz`, the pybedtools
  intervals, and the pseudobulk VCFs follow it, while everything internal stays
  chr-prefixed.

### Changed
- **Breaking**: `reference_version` is required on every sample-file record, and a run keeps
  only the records whose build matches the config's `reference_version`. One file may hold
  several builds; a run never mixes them. An unset config `reference_version` is now an
  error (it was a warning); an unrecognized value warns and is used verbatim.
- **Breaking**: the TSV sheet is a flat encoding of the JSON schema. Columns are the record
  keys and each input is its own `files.<key>` column, replacing `SAMPLE` / `REP_ID` /
  `PATH_to_bam` / `PATH_to_barcodes` / `PATH_to_10x_ranger`. TSV sheets can now carry remote
  URLs and name single-cell files explicitly.
- **Breaking**: `chromosomes` is validated against `genome_size` at DAG build; a chromosome
  the genome lacks is an error, not a silent drop. Every later step takes the list as
  validated.
- **Breaking**: `pseudobulk_snp_statistics.tsv` is no longer written, and its `report()`
  entry is gone from the single-cell DAG. `compute_snp_statistics` moved to
  `script_utils/annotate_snp_utils.py` and is marked legacy.
- The window build derives its contigs from `genome_size` instead of a hardcoded
  `chr1..chr22 + X/Y`, so a genome with more than 22 autosomes is no longer truncated.
- `parse_genetic_map` no longer guesses which numeric contig is X. A requested chromosome
  still missing after relabeling is an error naming it, where it was a warning; the
  numbering is consulted only when the map does not already label every chromosome.
- `do_repliseq` reads the canonical reference version, so an alias such as
  `reference_version: GRCh38` now enables Repli-seq RT correction; previously the alias was
  unrecognized and it stayed off.

### Fixed
- Contigs were renamed on disk for any reference without a `chr` prefix: `read_BED` and
  `read_VCF` normalize on ingest, and several scripts wrote the rewritten name back out.
  The prefix is now converted only where an external tool matches our output against the
  alignment or the reference.
- `build_window_bed` read the segment BED twice in two different naming styles, so on a
  bare-contig genome every window was dropped as off-segment.
- `rd_correct` (prebuilt `window_bed`, mosdepth output) and `atac_fragments_to_bb`
  (10x fragment files) now normalize contig names on read.
- `mappability_bed` intervals are renamed to the genome's convention before bedtools sees
  them, streamed per interval; previously a mismatch failed silently or with an opaque error.
- `gtf_file` contig names were passed through verbatim and joined against chr-normalized
  SNPs, so an Ensembl GTF annotated every SNP as `intergenic` with no warning and
  `gene_aware_binning` quietly degraded.
- A `files` entry set to `null` in a JSON sample file became the literal path `"None"` and
  failed DAG building; the key is now dropped, matching an empty TSV cell.
- A record missing a required key raised a bare `KeyError` from
  `validate_sample_file.py`, which subsetted by build before validating.

### Removed
- `CHR_STYLE_REFVERS`, `REFVER2SEXCHROM`, `get_standard_chroms`, `TSV_REQUIRED_COLUMNS`,
  and the TSV `PATH_to_10x_ranger` expansion.
- A duplicate `adaptive_dot_size` in `script_utils/utils.py`; every caller already used
  `cnplot`'s, which has the same signature.
- Three `chrom=` rule params that no shell referenced.

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
- Multi-replicate allele and depth consolidation onto one shared SNP/bin grid.

#### `bulk_genotyping`
- Bulk assays (WGS/WGS-lr/WES) on one shared window/bin grid, WES handled like WGS.
- bcftools genotyping from matched normals and het-SNP pileup; mosdepth depth with
  GC/mappability/replication-timing correction; adaptive binning to RDR/BAF matrices
  for HATCHet3.

#### `single_cell_genotyping`
- scRNA, scATAC, VISIUM/VISIUM3prime assays, including 10x Epi Multiome pairs.
- cellsnp-lite pseudobulk genotyping and per-cell pileup to allele and native-count
  matrices for CalicoST.

#### `copytyping_preprocess`
- Aggregate single-cell / spatial allele and native counts onto a pre-computed phased
  VCF and bin grid (no genotyping or phasing) for CalicoST.
