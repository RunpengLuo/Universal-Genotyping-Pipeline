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

### Changed
- **Breaking**: `reference_version` is required on every sample-file record, and a run keeps
  only the records whose build matches the config's `reference_version`. One file may hold
  several builds; a run never mixes them. An unset config `reference_version` is now an
  error (it was a warning); an unrecognized value warns and is used verbatim.
- **Breaking**: the TSV sheet is a flat encoding of the JSON schema. Columns are the record
  keys and each input is its own `files.<key>` column, replacing `SAMPLE` / `REP_ID` /
  `PATH_to_bam` / `PATH_to_barcodes` / `PATH_to_10x_ranger`. TSV sheets can now carry remote
  URLs and name single-cell files explicitly.

### Fixed
- A `files` entry set to `null` in a JSON sample file became the literal path `"None"` and
  failed DAG building; the key is now dropped, matching an empty TSV cell.
- A record missing a required key raised a bare `KeyError` from
  `validate_sample_file.py`, which subsetted by build before validating.

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
