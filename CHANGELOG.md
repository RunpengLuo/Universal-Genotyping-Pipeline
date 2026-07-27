# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-07-26

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
