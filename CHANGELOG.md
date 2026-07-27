# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-07-26

Initial pre-release.

### Added
- Snakemake pipeline for SNP genotyping, phasing, and allele counting, with three
  workflow modes: `bulk_genotyping`, `single_cell_genotyping`, `copytyping_preprocess`.
- Bulk assays (WGS/WGS-lr/WES) on one shared window/bin grid; WES handled identically
  to WGS.
- Single-cell / spatial assays (scRNA, scATAC, VISIUM, VISIUM3prime), including 10x Epi
  Multiome pairs.
- Genotyping via bcftools (bulk) and cellsnp-lite (single-cell pseudobulk), plus a
  tumor-only beta-binomial HMM caller when no matched normal is present.
- Phasing via Eagle2, SHAPEIT5, or LongPhase; `het_snp_vcf` short-circuits genotyping
  (and phasing when already phased).
- Bulk het-SNP read counting via `bcftools mpileup` aligned to the phased VCF.
- Read-depth counting (mosdepth) with LOWESS/median GC/MAP/replication-timing bias
  correction, adaptive binning, and RDR/BAF output matrices for HATCHet3 / Copy-typing
  / CalicoST.
- JSON sample file (with legacy TSV support) and remote (`http(s)`) inputs via
  Snakemake storage.
- Documentation: per-mode tutorials, config/output reference, sample-sheet spec, and
  external-resource catalog.
