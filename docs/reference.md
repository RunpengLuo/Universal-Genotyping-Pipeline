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
| `params_process_anndata` | `process_rna_anndata` | `gene_id_colname`, `min_frac_barcodes` |
| `params_phase_and_concat` | `phase_and_concat_{bulk,nonbulk}` | `min_depth` (bulk), `gamma` (bulk), `exon_only` |
| `params_mosdepth` | `run_mosdepth` | `read_quality`, `extra_params` |
| `params_count_reads` | `rd_correct` | `gc_correct`, `gc_correct_method` (`lowess`/`median`), `rt_correct`, `samplesize`, `routlier`, `doutlier`, `min_mappability` |
| `params_combine_counts` | `combine_counts`, `combine_counts_nonbulk` | `min_switchprob`, `nu`, `switchprob_ps`, `min_snp_reads`, `min_snp_per_block`, `gene_aware_binning`, `nsnp_multi` (sc only), `max_blocksize` (bulk only), `median_normalization` (bulk only), `rdr_outlier_quantile` (bulk only), `phase_flip_test` (bulk only), `phase_flip_epsilon` (bulk only), `phase_flip_alpha` (bulk only) |
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

### Allele Counts (`allele_dir/{stream_or_assay}/`)

Bulk writes ONE joint set under `allele_dir/{stream}/` (`stream` = `bulkWGS` for the WGS
family or `bulkWES`): all bulk assays segmented together; dense matrices, one pseudobulk
column per replicate. Single-cell writes per-assay sparse matrices under `allele_dir/{assay_type}/`.

- `snps.tsv.gz` — SNP annotations (chr, pos, ref/alt, phase, block).
- `snp.{Tallele,Aallele,Ballele}.npz` — allele count matrices (SNPs x samples/cells); dense for bulk, sparse for single-cell.
- `sample_ids.tsv` — sample metadata.
- `barcodes.tsv.gz` — cell barcode list (single-cell only).
- `barcodes.full.tsv.gz` — 2-col `REP_ID`/`BARCODE` mapping in matrix-column order (single-cell only).
- `unique_snp_ids.npy` — SNP identifiers as `{chr}_{pos}` (single-cell only).

### AnnData (`bb_dir/{assay_type}/`)

- `{assay_type}.h5ad` — AnnData with cells x features (single-cell only; produced by `process_anndata`).

### Final Bins — bulk: `bb_dir/{stream}/`; single-cell & copytyping: per-assay `bb_dir/{assay_type}/`

Common:
- `sample_ids.tsv` — sample metadata. For bulk and single-cell it carries an `assay_type` column and its row order matches the (sample / replicate×assay) matrix columns.

**Bulk (`bulk_genotyping`):** all bulk assays jointly segmented on one shared bin grid, written under `bb_dir/{stream}/` (`stream` = `bulkWGS` or `bulkWES`).
- `bb.tsv.gz` — bin annotations (one shared grid for all bulk assays).
- `bb.{Tallele,Aallele,Ballele,depth,rdr}.npz` — allele, depth, and RDR matrices. Columns concatenate all bulk samples (per assay, normal first); `rdr` holds the tumor columns only, normalized per assay against that assay's own normal. (BAF is not stored — derive it from `Ballele`/`Tallele`.)

**Single-cell (`single_cell_genotyping`):** all of the sample's non-bulk assays are jointly segmented on one shared grid (one pseudobulk column per replicate×assay). Everything lives under `bb_dir/{assay_type}/`; the shared grid and combined sample sheet are duplicated into each sub-dir.
- `{assay_type}/bb.tsv.gz` — the one shared bin grid (joint across assays; identical copy in each sub-dir).
- `{assay_type}/sample_ids.tsv` — one row per replicate×assay (identical copy in each sub-dir).
- `{assay_type}/bb.{Tallele,Aallele,Ballele}.npz` — per-assay allele count matrices (bins × cells) on the shared grid. (BAF is not stored — derive it from `Ballele`/`Tallele`.)
- `{assay_type}/bb.Xcount.npz` — per-assay native-count matrix per bb bin (bins × cells, sparse int32), same shape/column order as `{assay_type}/bb.{T,A,B}allele.npz`. **scATAC**: deduped ATAC fragment counts (each fragment counted once by its midpoint, from raw fragments). **scRNA / VISIUM / VISIUM3prime**: UMI counts, summed from the `process_rna_anndata` h5ad (each gene assigned to its largest-overlap bin, as in copytyping).
- `{assay_type}/multi_snp.tsv.gz`, `{assay_type}/multi_snp.{Tallele,Aallele,Ballele}.npz` — per-assay multi-SNP diagnostic groups.
- `{assay_type}/barcodes.tsv.gz`, `{assay_type}/barcodes.full.tsv.gz` — per-assay, copied from `allele_dir`.

**Copytyping (`copytyping_preprocess`):** per assay under `bb_dir/{assay_type}/`.
- `cnv_segments.tsv` — BB block annotations.
- `bb.{Xcount,Tallele,Aallele,Ballele}.npz` — per-block count matrices.
- `barcodes.tsv.gz`, `barcodes.full.tsv.gz` — copied from `allele_dir`.

### QC (`qc_dir/`, flat: `<stage>.<name>.<assay>.<run_id>.<ext>`)

All QC plots are written flat in `qc_dir`, with the pipeline stage, assay (or `bulk`), and run_id encoded in the filename:
- `phase_and_concat.snp_allele_freq.{assay}.{run_id}.pdf` — SNP allele frequency: page 1 = unphased REF/TOTAL, page 2 = phased B-allele frequency.
- `phase_and_concat.snp_depth_hist.{assay}.{run_id}.pdf` — SNP depth histogram.
- `rd_correction.rd_correct.{assay}.{run_id}.pdf` — read depth bias correction (bulk WGS/WES): before/after genome-wide RD scatter + GC/MAP/RT diagnostics.
- `combine_counts.combine_counts.bulk.{run_id}.pdf` — adaptive binning QC (bulk): segment-length & genes/segment histograms, per-sample count histograms, bin-level BAF and per-sample RDR/BAF scatter.
- `combine_counts.af_B_bb.pseudobulk.{assay}.{run_id}.pdf` (and `…af_B_multi-snp…`) — bin-level pseudobulk BAF (single-cell).
- `cnv_segmentation.af_cnv-B_{assay}.{assay}.{run_id}.pdf` — SNP and BB-level BAF (copytyping only).

---

## Key TSV Columns

### `snps.tsv.gz`

`#CHR`, `POS`, `POS0`, `START`, `END`, `GT`, `PHASE` (0 = B-allele is ALT, 1 = B-allele is REF), `region_id`, `feature_id` (`;`-joined list of all overlapping GTF genes, `intergenic` if none; GTF-derived for every assay), `feature_type` (exon/intron/intergenic). Bulk also carries `PS` (phase set) when the phaser (longphase) emits it.

### `bb.tsv.gz`

`#CHR`, `START`, `END`, `#SNPS`, `region_id`, `switchprobs`, `feature_id` (deduped `;`-joined union of the bin's SNP genes). `cnv_segments.tsv` carries the same `feature_id` column.

### `multi_snp.tsv.gz`

`#CHR`, `START`, `END`, `START0`, `END0`, `region_id`, `feature_id`, `feature_type`, `#SNPS`, `BLOCKSIZE`, `multi_id`, `switchprobs`.

### `sample_ids.tsv`

`SAMPLE` (`{patient_id}_{rep_id}`), `SAMPLE_NAME` (patient ID), `REP_ID`, `sample_type`.

### `barcodes.tsv.gz`

Single column, no header. Each row: `{BARCODE}_{REP_ID}`.

### `barcodes.full.tsv.gz`

Two columns with header: `REP_ID`, `BARCODE`. The `BARCODE` column carries the same suffixed `{BARCODE}_{REP_ID}` form as `barcodes.tsv.gz`, in matching row order. Used by QC plotters to group cells per replicate without parsing the suffix.
