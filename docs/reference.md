# Input / Output Reference

## Sample File

TSV with one row per replicate. Template at `resources/templates/samples.tsv`.

| Column | Required | Description |
|--------|----------|-------------|
| `SAMPLE` | Yes | Patient ID (matched by `sample_id` in config). |
| `REP_ID` | Yes | Unique per row; multiome pairs share a `REP_ID`. |
| `assay_type` | Yes | `bulkWGS`, `bulkWGS-lr`, `bulkWES`, `scATAC`, `scRNA`, `VISIUM`, or `VISIUM3prime`. |
| `sample_type` | Yes | `normal` or `tumor`. |
| `PATH_to_bam` | Yes | Path to `.bam` file. |
| `PATH_to_barcodes` | Non-bulk | Path to `barcodes.tsv.gz`. |
| `PATH_to_10x_ranger` | Non-bulk | Path to Cell Ranger / Space Ranger `outs/` directory. |

---

## Config Keys

Defaults live in `config/config.yaml`. A starting template for user runs is at `resources/templates/config.yaml`. Override individual keys via `--config key=value`.

### Top-level

| Key | Required | Description |
|-----|----------|-------------|
| `workflow_mode` | Yes | `bulk_genotyping` \| `single_cell_genotyping` \| `copytyping_preprocess`. |
| `assay_types` | Yes | List of assay types to run (e.g. `["bulkWGS"]`, `["scRNA","scATAC"]`). |
| `sample_id` | Yes | Selects which `SAMPLE` from the sample sheet to process. |
| `sample_file` | Yes | Path to `samples.tsv`. |
| `chromosomes` | Yes | List of chromosomes (default `[1..22]`; X handled separately by phasing tools). |
| `reference_version` | Yes | `hg19` \| `hg38` \| `chm13v2`. |
| `reference` | Yes | Genome FASTA. |
| `genome_size` | Yes | Two-column `chrom\tsize` file. |
| `region_bed` | Yes | Whitelist regions (e.g. autosomes minus blacklist). |
| `window_bed` | Bulk | Pre-built window BED with GC/MAP/REPLI columns. Build via `resources/scripts/build_{wgs,wes}_window_bed.py`. |
| `blacklist_bed` | Optional | ENCODE-style blacklist; pre-built at `resources/data/hg38-blacklist.v2.bed.gz`. |
| `gene_blacklist_file` | Optional | Genes to exclude from AnnData (single-cell). |
| `gtf_file` | Yes | Gene annotation GTF (gzipped). |
| `snp_panel` | Genotyping | Population SNP VCF. |
| `snp_targets` | Bulk genotyping | Per-chromosome position files; build via `resources/scripts/build_snp_targets.sh`. |
| `phasing_panel` | eagle/shapeit | Per-chromosome BCF reference panel directory. |
| `phaser` | Genotyping | `eagle` \| `shapeit` \| `longphase`. |
| `gmap_path` | eagle/shapeit | Full gmap path; use `{chrname}` placeholder for per-chrom files (SHAPEIT5) or a literal path for the single-file case (Eagle2). |
| `het_snp_vcf` | copytyping_preprocess | Pre-computed phased het SNP VCF (e.g. from a prior bulk run). |
| `bb_file` | copytyping_preprocess | Pre-computed BB block annotations TSV. |

### Parameter blocks

| Block | Used by | Keys |
|-------|---------|------|
| `params_cellsnp_lite` | `genotype_snps_pseudobulk_mode1b`, `pileup_snps_*` | `UMItag`, `cellTAG`, `minMAF_genotype`, `minCOUNT_genotype`, `minMAF_pileup`, `minCOUNT_pileup` |
| `params_bcftools` | `genotype_snps_bulk` | `min_mapq`, `min_baseq`, `min_dp`, `max_depth`, `min_qual` |
| `params_annotate_snps` | `annotate_snps_pseudobulk` | `min_het_reads`, `min_hom_dp`, `min_vaf_thres`, `filter_nz_OTH`, `filter_hom_ALT` |
| `params_longphase` | `phase_snps_longphase` | `min_mapq`, `extra_params` (`--pb` or `--ont`) |
| `params_process_anndata` | `process_rna_anndata`, `process_atac_fragments` | `gene_id_colname`, `min_frac_barcodes`, `tile_width` |
| `params_phase_and_concat` | `phase_and_concat_{bulk,single_cell}` | `min_depth`, `gamma`, `exon_only` |
| `params_mosdepth` | `run_mosdepth` | `read_quality`, `extra_params` |
| `params_count_reads` | `rd_correct` | `gc_correct`, `gc_correct_method` (`lowess`/`median`), `rt_correct`, `samplesize`, `routlier`, `doutlier`, `min_mappability` |
| `params_combine_counts` | `combine_counts`, `combine_counts_nonbulk` | `min_switchprob`, `nu`, `switchprob_ps`, `nsnp_multi` (sc only), `min_snp_reads`, `min_snp_per_block`, `max_blocksize` (bulk only), `median_normalization`, `rdr_outlier_quantile` (bulk only), `phase_flip_test`, `phase_flip_epsilon`, `phase_flip_alpha` |
| `threads` | All multi-thread rules | `genotype`, `phase`, `pileup`, `mosdepth` |

### Output directories

`snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `bb_dir`, `qc_dir`, `log_dir` — all relative to `snakemake --directory`.

---

## Outputs

All output directories (`snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `bb_dir`, `qc_dir`, `log_dir`) are set in `config.yaml`.

### SNP Genotyping (`snp_dir/`)

- `chr{chrname}.vcf.gz` — per-chromosome VCF of bi-allelic SNPs.
- `pseudobulk_{modality}/cellSNP.*` — pseudobulk pileup (single-cell only).

### Phasing (`phase_dir/`)

- `chr{chrname}.vcf.gz` — per-chromosome phased VCF.
- `phased_het_snps.vcf.gz(.tbi)` — genome-wide phased het SNPs.
- `germline_snp_statistics.tsv` — per-chromosome counts of het_phased, het_unphased, hom_alt, hom_ref SNPs.
- `genetic_map.tsv.gz` — cM positions per SNP (eagle/shapeit only).

### Pileup (`pileup_dir/`)

One subdirectory per `{assay_type}_{rep_id}` with cellsnp-lite output.

### Allele Counts (`allele_dir/{assay_type}/`)

- `snps.tsv.gz` — SNP annotations (chr, pos, ref/alt, phase, block).
- `snp.{Tallele,Aallele,Ballele}.npz` — sparse allele count matrices (SNPs x samples/cells).
- `sample_ids.tsv` — sample metadata.
- `barcodes.tsv.gz` — cell barcode list (single-cell only).
- `barcodes.full.tsv.gz` — 2-col `REP_ID`/`BARCODE` mapping in matrix-column order (single-cell only).
- `unique_snp_ids.npy` — SNP identifiers as `{chr}_{pos}` (single-cell only).
- `cell_snp_Aallele.npz`, `cell_snp_Ballele.npz` — per-cell allele count matrices (single-cell only).

### AnnData (`bb_dir/{assay_type}/`)

- `{assay_type}.h5ad` — AnnData with cells x features (single-cell only; produced by `process_anndata`).

### Final Bins (bulk: `bb_dir/bulkWGS/`, or `bb_dir/bulkWES/` for a WES-only run; single-cell/copytyping: `bb_dir/{assay_type}/`)

Common outputs across all modes:
- `sample_ids.tsv` — sample metadata. For bulk this carries an `assay_type` column and its row order matches the matrix columns.

**Bulk (`bulk_genotyping`):** all bulk assays are jointly segmented on one shared bin grid under `bb_dir/bulkWGS/` (WGS family) or `bb_dir/bulkWES/` (WES-only).
- `bb.tsv.gz` — bin annotations (one shared grid for all bulk assays).
- `bb.{Tallele,Aallele,Ballele,baf,depth,rdr}.npz` — allele, depth, and RDR matrices. Columns concatenate all bulk samples (per assay, normal first); `rdr` holds the tumor columns only, normalized per assay against that assay's own normal.

**Single-cell (`single_cell_genotyping`):**
- `bb.tsv.gz` — bin annotations.
- `bb.{Tallele,Aallele,Ballele}.npz` — allele count matrices.
- `multi_snp.tsv.gz` — multi-SNP group annotations.
- `multi_snp.{Tallele,Aallele,Ballele}.npz` — multi-SNP allele count matrices.
- `barcodes.tsv.gz`, `barcodes.full.tsv.gz` — copied from `allele_dir`.

**Copytyping (`copytyping_preprocess`):**
- `cnv_segments.tsv` — BB block annotations.
- `bb.{Xcount,Tallele,Aallele,Ballele}.npz` — per-block count matrices.
- `barcodes.tsv.gz`, `barcodes.full.tsv.gz` — copied from `allele_dir`.

### QC (bulk: `qc_dir/bulkWGS/` or `qc_dir/bulkWES/`; single-cell/copytyping: `qc_dir/{assay_type}/`)

Plots from bias correction, phasing, and binning steps:
- `rd_correct.{run_id}.pdf` — read depth bias correction (bulk WGS/WES). One page per sample with before/after correction genome-wide RD scatter, followed by GC correction diagnostic pages.
- `snp_allele_freq.{run_id}.pdf` — SNP-level allele frequency: page 1 = unphased REF/TOTAL, page 2 = phased B-allele frequency.
- `snp_depth_hist.{run_id}.pdf` — SNP depth histogram.
- `combine_counts.{run_id}.pdf` — adaptive binning QC: bin-level BAF + per-sample RDR/BAF (bulk).
- `af_cnv-B_{assay_type}.{run_id}.pdf` — SNP and BB-level BAF (copytyping only).

---

## Key TSV Columns

### `snps.tsv.gz`

`#CHR`, `POS`, `POS0`, `START`, `END`, `GT`, `PHASE` (0 = B-allele is ALT, 1 = B-allele is REF), `region_id`, `feature_id`, `feature_type` (exon/intron/intergenic).

### `bb.tsv.gz`

`#CHR`, `START`, `END`, `#SNPS`, `region_id`, `switchprobs`.

### `multi_snp.tsv.gz`

`#CHR`, `START`, `END`, `START0`, `END0`, `region_id`, `feature_id`, `feature_type`, `#SNPS`, `BLOCKSIZE`, `multi_id`, `switchprobs`.

### `sample_ids.tsv`

`SAMPLE` (`{patient_id}_{rep_id}`), `SAMPLE_NAME` (patient ID), `REP_ID`, `sample_type`.

### `barcodes.tsv.gz`

Single column, no header. Each row: `{BARCODE}_{REP_ID}`.

### `barcodes.full.tsv.gz`

Two columns with header: `REP_ID`, `BARCODE`. The `BARCODE` column carries the same suffixed `{BARCODE}_{REP_ID}` form as `barcodes.tsv.gz`, in matching row order. Used by QC plotters to group cells per replicate without parsing the suffix.
