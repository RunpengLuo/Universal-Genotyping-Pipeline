# Reference

## Table of Contents
- [Environments](#environments)
- [Configuration](#configuration)
  - [Input Data](#input-data)
  - [Parameters](#parameters)
- [Outputs](#outputs)
  - [Genomic unit levels](#genomic-unit-levels)
  - [Final bins](#final-bins)
  - [Intermediates](#intermediates)
  - [TSV columns](#tsv-columns)
  - [QC](#qc-qc_dir)

---

## Environments

Pipeline dependencies live under `workflow/envs/`:

| Environment | Purpose |
|-------------|---------|
| `base.yaml` | Python scientific stack (used by every mode). |
| `bcftools.yaml` | bcftools/tabix: bulk genotyping + het-SNP pileup. |
| `eagle.yaml` | Eagle2 phasing. |
| `shapeit.yaml` | SHAPEIT5 phasing. |
| `longphase.yaml` | LongPhase (long-read) phasing. |
| `cellsnp.yaml` | cellsnp-lite: single-cell genotyping + pileup. |
| `mosdepth.yaml` | mosdepth read-depth counting (bulk). |
| `ucsc.yaml` | UCSC tools (`bigWigToBedGraph`, `liftOver`) for the Repli-seq track. |

---

## Configuration
Defaults in `config/config.yaml`, template in [templates](../resources/templates/). Override with `--config key=value`. Refer to **[spec](sample_sheet.md)** and [templates](../resources/templates/) for sample sheet.

### Input Data

| Key | Required | Description |
|-----|----------|-------------|
| `workflow_mode` | Yes | `bulk_genotyping` \| `single_cell_genotyping` \| `copytyping_preprocess`. |
| `assay_types` | Yes | Assay types to run, e.g. `["bulkWGS"]`, `["scRNA","scATAC"]`. |
| `sample_id` | Yes | Which `sample_id` of the sample file to process. |
| `sample_file` | Yes | Path to `samples.json`. |
| `chromosomes` | Yes | Chromosomes to run; default `[1..22]`. |
| `remote_mode` | Optional | Remote input handling: `storage` (default; download whole file via Snakemake storage) or `stream` (read URLs directly, fetching only `chromosomes`). |
| `species` | Yes | `human` (default) or `mouse`. |
| `reference_version` | Yes | Reference version to select samples. See [Reference version](sample_sheet.md#reference-version). |
| `reference` | Yes | Genome FASTA. |
| `genome_size` | Yes | Two-column `chrom\tsize` genome size file. |
| `region_bed` | Yes | BED file listing whitelist chromosome arms. Column 4 (BED NAME), when present, is the `region_id`; a BED3 gets `CHR:START-END`. |
| `segment_bed` | Optional | BED file listing genomic segments separated by novel adjacency. Column 4 (BED NAME), when present, is the `seg_id`; a BED3 gets `CHR:START-END`. Unset -> `region_bed`, one segment per arm. |
| `window_bed` | Optional | Pre-built window BED with read depth covariates. [pre-built](`resources/data/windows.1kbp.{hg19,hg38,chm13v2}.bed.gz`). |
| `gtf_file` | Yes | Gene annotation GTF (gzipped). |
| `mappability_bed` | Optional | BED mappability track (4th column = score). |
| `blacklist_bed` | Optional | ENCODE-style blacklist; pre-built at `resources/data/hg38-blacklist.v2.bed.gz`. |
| `gene_blacklist_file` | Optional | Genes to exclude from AnnData (single-cell). |
| `snp_panel` | Genotyping | Population SNP VCF, bgzipped and indexed (`.vcf.gz` + `.tbi`/`.csi`). Bulk passes it to `bcftools mpileup -T` (positions only, panel alleles ignored); single-cell to `cellsnp-lite -R`. |
| `phaser` | Genotyping | `eagle` \| `shapeit` \| `longphase`. |
| `phasing_panel` | eagle/shapeit | Per-chromosome BCF reference panel directory. |
| `gmap_path` | eagle/shapeit | Genetic map; `{chrname}` placeholder for per-chromosome maps (SHAPEIT5), literal path for a single map (Eagle2). |
| `genotype_dataset_ids` | Optional | `dataset_id`s piled up to call germline SNPs. Empty -> auto (normal before tumor, short-read before long-read). >1 are pooled in one `mpileup` and must share an `@RG SM` tag. |
| `phase_dataset_ids` | Optional | datasets used for long-read phasing inputs |
| `het_snp_vcf` | Optional; required for `copytyping_preprocess` | Pre-computed gHET VCF. |
| `het_snp_vcf_phased` | Optional | Default `true`: the input `het_snp_vcf` is phased or not. |
| `bb_file` | copytyping_preprocess | Pre-computed bb annotations TSV. |

> [!IMPORTANT]
> - `genome_size` defines the contig naming convention in `reference` and input alignment files.
> - `snp_panel`, `phasing_panel` and `het_snp_vcf` must follow the same naming
> convention as `genome_size`.
> - Final outputs always use chr-prefix contig naming regardless of input convention.

### Parameters
#### `params_build_windows`
Used by the window-BED build (`build_windows.smk`), which runs in **every** mode.
`build_segment_bed` first stamps each `segment_bed` segment with the arm it sits in and
subtracts the blacklist into `aux/segment.bed` (region_id + seg_id). One window BED is
then tiled off it, per segment row, so no window and no bin spans a segment bound:
`aux/windows.bed.gz`, shared by every assay of the run.

The windows are the fixed bins of `build_adaptive_bins` in both modes. Bulk additionally
counts them with mosdepth and bias-corrects them in `rd_correct`, which is the only
consumer of the GC/MAP/REPLI columns; a single-cell window BED carries just
`#CHR START END region_id seg_id`, and the Repli-seq fetch is skipped (as it is for any run with `rt_correct: false`).

| Field | Description |
|---|---|
| `window_size` | Window size (bp); fixed tiling of the segment BED, shared by every assay of the run. |

#### `params_bcftools`
Used by `genotype_snps_bulk` and `pileup_snps_bulk_bcftools_chrom` (bulk het-SNP read
counting). Both call one chromosome per job and restrict it with `--regions`.

| Field | Description |
|---|---|
| `min_mapq` | Skip alignments below this mapping quality (genotype + pileup). |
| `min_baseq` | Skip bases below this base quality (genotype + pileup). |
| `min_dp` | Minimum depth to keep a site (genotype only). |
| `max_depth` | Per-file depth cap in `mpileup` (genotype + pileup). |
| `min_qual` | Minimum variant QUAL (genotype only). |
| `extra_params` | Extra `mpileup` flags applied to bulk genotyping + pileup; use a read-type preset, e.g. `-X ont-sup` (ONT), `-X pacbio-ccs`, or `""`/`-X illumina` for short-read. |

#### `params_cellsnp_lite`
Used by `genotype_snps_pseudobulk_mode1b`, `pileup_snps_*` (single-cell).

| Field | Description |
|---|---|
| `UMItag` | UMI tag: `Auto` \| `None` \| a BAM tag. |
| `cellTAG` | Cell-barcode tag: `CB` \| `None` \| a BAM tag. |
| `minMAF_genotype` | Minimum minor-allele frequency when genotyping. |
| `minCOUNT_genotype` | Minimum aggregate count when genotyping. |
| `minMAF_pileup` | Minimum minor-allele frequency when piling up. |
| `minCOUNT_pileup` | Minimum aggregate count when piling up. |

#### `params_annotate_snps`
Used by `genotype_snps_no_normal` (single-cell), which genotypes from read counts because
no matched normal exists to call against.

| Field | Description |
|---|---|
| `min_het_reads` | Minimum reads on *each* allele for a het call. |
| `min_hom_dp` | Minimum depth for a hom call. |
| `min_vaf_thres` | Het VAF must lie in `[min_vaf_thres, 1 - min_vaf_thres]`. |
| `filter_nz_OTH` | Drop SNPs with non-zero OTH (non-ref, non-alt) counts. |
| `filter_hom_ALT` | Drop hom-ALT SNPs. |

#### `params_longphase`
Used by `phase_snps_longphase`.

| Field | Description |
|---|---|
| `min_mapq` | Skip alignments below this mapping quality. |
| `extra_params` | Read chemistry: `--pb` (PacBio) \| `--ont` (ONT). |

#### `params_process_anndata`
Used by `process_rna_anndata` (single-cell RNA / spatial).

| Field | Description |
|---|---|
| `gene_id_colname` | Gene-id column joined against the GTF. |
| `min_frac_barcodes` | Drop a gene expressed in fewer than this fraction of barcodes. |

#### `params_phase_and_concat`
Used by `phase_and_concat_{bulk,nonbulk}`.

| Field | Description |
|---|---|
| `min_depth` | Minimum depth in every sample to keep a SNP (bulk). |
| `gamma` | Credible-interval level of the balanced-het test on the normal; a SNP is kept when its beta posterior interval covers 0.5 (bulk). |
| `exon_only` | Keep exonic SNPs only. |

#### `params_mosdepth`
Used by `run_mosdepth` (bulk).

| Field | Description |
|---|---|
| `read_quality` | Skip alignments below this mapping quality. |
| `extra_params` | Extra mosdepth flags. |

#### `params_count_reads`
Used by `rd_correct` (bulk); HMMcopy-style bias correction.

| Field | Description |
|---|---|
| `rd_correct_method` | `lowess` \| `median`. |
| `gc_correct` | Model GC content as a covariate during fitting. |
| `rt_correct` | Model replication-timing as a covariate during fitting. |
| `samplesize` | Number of sampled bins during fitting. |
| `routlier` | Upper quantile of read count dropped as outlier. |
| `doutlier` | Top/bottom quantile of the GC/mappability domain dropped as outlier. |
| `min_mappability` | Drop bins with mappability below this cutoff. |

#### `params_combine_counts`
Used by `combine_counts` (bulk) and `combine_counts_nonbulk` (single-cell).

| Field | Description |
|---|---|
| `min_snp_reads` | Reads every tumor column needs to close a bb; a list sweeps `MSR{msr}/`. |
| `min_snp_per_bin` | SNPs needed to close a bb. |
| `gene_aware_binning` | Grow bbs by whole genes; never cut inside one. |
| `nu` | Haldane scale turning cM distance into a switch probability. |
| `min_switchprob` | Floor on that switch probability. |
| `switchprob_ps` | Switch probability within one phase set (`PS`); ~0.5 across sets. |
| `nsnp_multi` | SNPs per multi-SNP diagnostic group. |
| `max_blocksize` | Span (bp) past which `min_snp_per_bin` is waived; `0` = no cap (bulk). |
| `rdr_normalization` | Bulk RDR denominator: `auto` (base else median), `median`, `normal` (base required). |
| `phase_flip_test` | Split a phase cluster failing the haplotype-flip test (bulk). |
| `phase_flip_epsilon` | Effect size of that test (bulk). |
| `phase_flip_alpha` | Significance level of that test (bulk). |

> [!NOTE]
> Adaptive binning merges consecutive windows, and closes a bb only when:
> - every tumor column has `min_snp_reads` SNP reads.
> - the bb holds `min_snp_per_bin` SNPs, or `[bulk]` spans `max_blocksize`.
> - the next window starts a new gene, under `gene_aware_binning`.
> - a bb never spans two `region_id`, `seg_id` or `PS` clusters.
> - `[bulk]` nor two `phase_cluster` clusters, under `phase_flip_test`.
> - a trailing run below threshold merges into the previous bb.
>
> bbs are then post-filtered:
> - `[bulk]` drop a bb whose BAF is NaN: no SNP, or a column with no read over them.
> - `[bulk]` drop a bb whose depth is NaN: every window below `min_mappability`, or
>   (`lowess`) outside the fitted range; or a multi-SNP group overlapping no window.
> - `[bulk]` drop a bb whose RDR is NaN: any depth NaN, a zero or NaN base column, or an
>   all-NaN tumor column.
> - `[nonbulk]` no filter; SNP-free bbs are kept, with all-zero allele rows and their
>   `Xcount`.

#### `threads`
Used by all multi-thread rules.

| Field | Description |
|---|---|
| `genotype` | Threads for genotyping; `bcftools mpileup` is single-threaded, so these size the `call` and `view` output compressors only. |
| `phase` | Threads for phasing. |
| `pileup` | Threads for the pileup step (bulk: the `bgzip` writing the counts; single-cell: cellsnp-lite, which is genuinely parallel). |
| `mosdepth` | Threads for mosdepth. |

> [!NOTE]
> `bcftools` applies `--threads` to the output handle only, never to the BAM/CRAM
> readers, so a bulk `mpileup` runs on one core whatever this is set to. Bulk genotyping
> and pileup both fan out per chromosome instead, and small values here leave cores free
> for more concurrent chromosome jobs.

---

## Outputs

Set in `config.yaml`, relative to `snakemake --directory`:

| Key | Path | Contents |
|---|---|---|
| `snp_dir` | `snps` | Genotyped SNP VCFs. |
| `phase_dir` | `phase` | Phased VCFs, parsed genetic map. |
| `pileup_dir` | `pileup` | Per-dataset allele counts and read depth. |
| `allele_dir` | `allele` | The SNP grid and its allele matrices. |
| `bb_dir` | `bb` | The bbs and their matrices. |
| `qc_dir` | `qc` | One multi-page PDF per rule. |
| `log_dir` | `logs` | `{rule}/...`, one log per job. |
| `aux_dir` | `aux` | Segment BED, window BED, Repli-seq tracks. |
| `bench_dir` | `benchmarks` | `{rule}/...`, runtime and `max_rss` per job. |

### Genomic unit levels

| Level | id | Description |
|---|---|---|
| Region | `region_id` | Chromosome arm, from `region_bed`. Carried for RDR and QC. |
| Segment | `seg_id` | From `segment_bed`, the hard bb bound. Blacklist pieces keep it. |
| Window | `bin_id` | `window_size` tile, shared by every assay. Never written out. |
| bb | `bb_id` | Merged windows; the feature axis of `bb.tsv.gz` and `bb.*.npz`. |

### Final bins

A `min_snp_reads` list writes one `MSR{msr}/` per value.

| Mode | Location |
|---|---|
| `bulk_genotyping` | `bb_dir/MSR{msr}/bulk/` |
| `single_cell_genotyping` | `bb_dir/MSR{msr}/{assay_type}/` |
| `copytyping_preprocess` | `bb_dir/{assay_type}/` |

**`bulk_genotyping`** - one bb set over every bulk assay; columns are samples.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations, the matrix row axis. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `bb.{depth,rdr}.npz` | Depth, and RDR for the tumor columns only. |
| `sample_ids.tsv` | One row per matrix column. |
| `multi_snp/bulk/` | Multi-SNP diagnostic groups, bb schema, outside `MSR{msr}/`. |

**`single_cell_genotyping`** - one bb set over every assay, copied into each subdir;
columns are cells.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations, the matrix row axis. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `bb.Xcount.npz` | Native counts: scATAC fragment midpoints, scRNA/VISIUM h5ad UMIs. |
| `sample_ids.tsv` | One row per dataset, not column-aligned. |
| `barcodes.tsv.gz` | The matrix column axis; each barcode is suffixed `_{dataset_id}_{assay_type}`. |
| `multi_snp/{assay_type}/` | Multi-SNP diagnostic groups, bb schema, outside `MSR{msr}/`. |

**`copytyping_preprocess`** - the given bb set, per assay; columns are cells.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations, the matrix row axis; the given `bb_file` plus `#SNPS`, `feature_id`, (RNA) `#feature`. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts, in matrix-column order. |
| `bb.Xcount.npz` | Native counts: scATAC fragment midpoints, scRNA/VISIUM h5ad UMIs. |
| `sample_ids.tsv` | One row per dataset, not column-aligned. |
| `barcodes.tsv.gz` | The matrix column axis; each barcode is suffixed `_{dataset_id}_{assay_type}`. |

### Intermediates

| Path | Contents |
|---|---|
| `snp_dir/chr{chrname}.vcf.gz` | Bi-allelic SNPs. |
| `snp_dir/pseudobulk_{modality}/` | cellsnp-lite pseudobulk output. |
| `phase_dir/chr{chrname}.vcf.gz` | Phased SNPs, concatenated to `phased_het_snps.vcf.gz`. |
| `phase_dir/genetic_map.tsv.gz` | Parsed genetic map (eagle/shapeit). |
| `pileup_dir/{assay_type}_{dataset_id}/` | Bulk `bcftools.counts.tsv.gz`, concatenated from per-chromosome `bcftools.counts.chr{chrname}.tsv.gz` (temporary); single-cell `cellSNP.*`. |
| `pileup_dir/{assay_type}/out_mosdepth/` | Per-dataset mosdepth (bulk). |
| `pileup_dir/bulk/window.dp.npz` | Corrected depth, windows x every bulk dataset. |
| `allele_dir/` | `snps.tsv.gz`, `snp.{T,A,B}allele.npz`, `sample_ids.tsv`, `barcodes.tsv.gz`. |
| `bb_dir/{assay_type}.h5ad` | Gene x cell AnnData (scRNA / spatial). |
| `aux_dir/segment.bed` | Arm-stamped segments, blacklist subtracted. |
| `aux_dir/windows.bed.gz` | The shared window BED. |
| `aux_dir/repliseq/` | Repli-seq tracks, lifted from hg19 when the run is not hg19. |

### TSV columns

| File | Columns |
|---|---|
| `snps.tsv.gz` | `#CHR POS POS0 START END GT PHASE region_id seg_id feature_id feature_type`; bulk adds `PS`. |
| `bb.tsv.gz` | `#CHR START END #SNPS region_id switchprobs feature_id`. |
| `sample_ids.tsv` | `SAMPLE sample_id dataset_id sample_type assay_type`; bulk adds `rdr_base_dataset_id`. |
| `germline_snp_statistics.tsv` | Per chromosome: het_phased, het_unphased, hom_alt, hom_ref. |
| `depth_statistics.tsv` | Per-dataset depth summary, every bulk dataset in one table. |

> [!NOTE]
> - `PHASE`: 0 = the B-allele is ALT, 1 = the B-allele is REF.
> - `feature_id`: `;`-joined overlapping GTF genes, `intergenic` if none.
> - `sample_ids.tsv` is column-aligned in bulk only; single-cell columns are cells.
> - `SAMPLE` is `{sample_id}_{dataset_id}`, plus `_{assay_type}` for a multiome pair.

### QC (`qc_dir/`)

One multi-page PDF per rule, flat:

| File | Contents |
|---|---|
| `phase_and_concat.{bulk_or_assay}.pdf` | SNP allele frequency and depth. |
| `rd_correction.bulk.pdf` | Depth before/after correction, GC/MAP/RT diagnostics. |
| `combine_counts.{bulk_or_assay}.MSR{msr}.pdf` | Binning QC, one per `min_snp_reads`. |
| `combine_counts_fixed_bins.{assay_type}.pdf` | SNP- and bb-level BAF. |
