# Input / Output Reference

## Sample File
Refer to spec **[sample_sheet.md](sample_sheet.md)** and [templates](../resources/templates/).

---

## Config Keys
Defaults in `config/config.yaml`, template in [templates](../resources/templates/). Override with `--config key=value`.

### Input Data

| Key | Required | Description |
|-----|----------|-------------|
| `workflow_mode` | Yes | `bulk_genotyping` \| `single_cell_genotyping` \| `copytyping_preprocess`. |
| `assay_types` | Yes | Assay types to run, e.g. `["bulkWGS"]`, `["scRNA","scATAC"]`. |
| `sample_id` | Yes | Which `sample_id` of the sample file to process. |
| `sample_file` | Yes | Path to `samples.json`. |
| `chromosomes` | Yes | Chromosomes to run; default `[1..22]`. |
| `reference_version` | Yes | `hg19` \| `hg38` \| `chm13v2` \| `mm10`. Another value runs, with a warning. |
| `reference` | Yes | Genome FASTA. |
| `genome_size` | Yes | Two-column `chrom\tsize` file. |
| `region_bed` | Yes | Whitelist regions. |
| `gtf_file` | Yes | Gene annotation GTF (gzipped). |
| `window_bed` | Bulk | Window BED with GC/MAP/REPLI columns; build via `resources/scripts/build_{wgs,wes}_window_bed.py`. |
| `blacklist_bed` | Optional | ENCODE-style blacklist; pre-built at `resources/data/hg38-blacklist.v2.bed.gz`. |
| `gene_blacklist_file` | Optional | Genes to exclude from AnnData (single-cell). |
| `snp_panel` | Genotyping | Population SNP VCF. |
| `snp_targets` | Bulk genotyping | Per-chromosome position files; build via `resources/scripts/build_snp_targets.sh`. |
| `phaser` | Genotyping | `eagle` \| `shapeit` \| `longphase`. |
| `phasing_panel` | eagle/shapeit | Per-chromosome BCF reference panel directory. |
| `gmap_path` | eagle/shapeit | Genetic map; `{chrname}` placeholder for per-chromosome maps (SHAPEIT5), literal path for a single map (Eagle2). |
| `genotype_dataset_ids` | Optional | `dataset_id`s piled up to call germline SNPs. Empty -> auto (normal before tumor, short-read before long-read). >1 are pooled in one `mpileup` and must share an `@RG SM` tag. |
| `phase_dataset_ids` | Optional | The one `dataset_id` `longphase` reads. Empty -> auto (long-read normal, else long-read tumor). Ignored by panel phasers. |
| `het_snp_vcf` | Optional; required for `copytyping_preprocess` | Pre-computed het SNP VCF. Set in any mode to skip genotyping. |
| `het_snp_vcf_phased` | Optional | Default `true`: the VCF is taken as phased, so phasing is skipped too. `false` phases it. Read only with `het_snp_vcf`. |
| `bb_file` | copytyping_preprocess | Pre-computed BB block annotations TSV. |

### Parameters

Defaults are those in `config/config.yaml`.

#### `params_cellsnp_lite`
Used by `genotype_snps_pseudobulk_mode1b`, `pileup_snps_*` (single-cell).

| Field | Default | Description |
|---|---|---|
| `UMItag` | `Auto` | UMI tag: `Auto` \| `None` \| a BAM tag. |
| `cellTAG` | `CB` | Cell-barcode tag: `CB` \| `None` \| a BAM tag. |
| `minMAF_genotype` | `0` | Minimum minor-allele frequency when genotyping. |
| `minCOUNT_genotype` | `2` | Minimum aggregate count when genotyping. |
| `minMAF_pileup` | `0` | Minimum minor-allele frequency when piling up. |
| `minCOUNT_pileup` | `1` | Minimum aggregate count when piling up. |

#### `params_bcftools`
Used by `genotype_snps_bulk`.

| Field | Default | Description |
|---|---|---|
| `min_mapq` | `20` | Skip alignments below this mapping quality. |
| `min_baseq` | `20` | Skip bases below this base quality. |
| `min_dp` | `5` | Minimum depth to keep a site. |
| `max_depth` | `1000` | Per-file depth cap in `mpileup`. |
| `min_qual` | `30` | Minimum variant QUAL. |

#### `params_annotate_snps`
Used by `annotate_snps_pseudobulk` (single-cell).

| Field | Default | Description |
|---|---|---|
| `min_het_reads` | `2` | Minimum reads on *each* allele for a het call. |
| `min_hom_dp` | `10` | Minimum depth for a hom call. |
| `min_vaf_thres` | `0.1` | Het VAF must lie in `[min_vaf_thres, 1 - min_vaf_thres]`. |
| `filter_nz_OTH` | `false` | Drop SNPs with non-zero OTH (non-ref, non-alt) counts. |
| `filter_hom_ALT` | `false` | Drop hom-ALT SNPs. |

#### `params_longphase`
Used by `phase_snps_longphase`.

| Field | Default | Description |
|---|---|---|
| `min_mapq` | `20` | Skip alignments below this mapping quality. |
| `extra_params` | `--pb` | Read chemistry: `--pb` (PacBio) \| `--ont` (ONT). |

#### `params_process_anndata`
Used by `process_rna_anndata` (single-cell RNA / spatial).

| Field | Default | Description |
|---|---|---|
| `gene_id_colname` | `gene_ids` | Gene-id column joined against the GTF. |
| `min_frac_barcodes` | `5e-3` | Drop a gene expressed in fewer than this fraction of barcodes. |

#### `params_phase_and_concat`
Used by `phase_and_concat_{bulk,nonbulk}`.

| Field | Default | Description |
|---|---|---|
| `min_depth` | `1` | Minimum depth in every sample to keep a SNP (bulk). |
| `gamma` | `0.05` | Credible-interval level of the balanced-het test on the normal; a SNP is kept when its beta posterior interval covers 0.5 (bulk). |
| `exon_only` | `false` | Keep exonic SNPs only. |

#### `params_mosdepth`
Used by `run_mosdepth` (bulk).

| Field | Default | Description |
|---|---|---|
| `read_quality` | `11` | Skip alignments below this mapping quality. |
| `extra_params` | `--no-per-base --fast-mode` | Extra mosdepth flags. |

#### `params_count_reads`
Used by `rd_correct` (bulk); HMMcopy-style bias correction.

| Field | Default | Description |
|---|---|---|
| `gc_correct` | `true` | Correct GC bias. |
| `gc_correct_method` | `median` | `lowess` \| `median`. |
| `rt_correct` | `true` | Also correct replication-timing bias (needs a `REPLI` column in `window_bed`). |
| `samplesize` | `50000` | Max ideal windows used to fit the correction curve. |
| `routlier` | `0.01` | Upper quantile of read count dropped as outlier. |
| `doutlier` | `0.001` | Top/bottom quantile of the GC/mappability domain dropped as outlier. |
| `min_mappability` | `0.9` | Minimum mappability for a window to be "ideal". |

#### `params_combine_counts`
Used by `combine_counts` (bulk) and `combine_counts_nonbulk` (single-cell).

| Field | Default | Description |
|---|---|---|
| `min_snp_reads` | `[500, 1000]` | SNP-read target per bin. A list sweeps binning: one job preprocesses once and writes one `MSR{msr}/` subdir per value. |
| `min_snp_per_block` | `2` | Minimum SNPs per bin. |
| `gene_aware_binning` | `true` | `true` = a bin grows by whole genes (GTF `feature_id`) and never splits one; `false` = window/SNP-granular binning. |
| `nu` | `1` | Scale of the Haldane map `(1 - exp(-2*nu*d)) / 2` turning cM distance into a switch probability. |
| `min_switchprob` | `1e-6` | Floor on that switch probability. |
| `switchprob_ps` | `5e-2` | Switch probability within one phase set (`PS`); across phase sets it is ~0.5. |
| `nsnp_multi` | `2` | SNPs per multi-SNP diagnostic group (single-cell). |
| `max_blocksize` | `500000` | Force a bin cut past this span, in bp (bulk; `0` = no cap). |
| `rdr_normalization` | `auto` | RDR denominator of a bulk tumor: `auto` = its `rdr_base_dataset_id` when set, else median; `median` = always median; `normal` = always `rdr_base_dataset_id`, error if any tumor lacks one. Resolved at parse time; a tumor's base should be a same-platform normal so RDR cancels platform bias. |
| `rdr_outlier_quantile` | `0.002` | Bins above the `1 - q` RDR quantile are set to NaN (bulk; `0` = off). |
| `phase_flip_test` | `true` | Test each phase-set block for a haplotype flip and split it (bulk). |
| `phase_flip_epsilon` | `0.05` | Effect size of that test (bulk). |
| `phase_flip_alpha` | `0.05` | Significance level of that test (bulk). |

#### `threads`
Used by all multi-thread rules.

| Field | Default | Description |
|---|---|---|
| `genotype` | `4` | Threads for genotyping. |
| `phase` | `4` | Threads for phasing. |
| `pileup` | `8` | Threads for cellsnp-lite pileup. |
| `mosdepth` | `4` | Threads for mosdepth. |

---

## Outputs

Directories (`snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `bb_dir`, `qc_dir`, `log_dir`) are set in `config.yaml`, relative to `snakemake --directory`. Which rule writes what, per mode: [workflow.md](workflow.md).

`.npz` are matrices: rows = SNPs or bins, columns = samples or cells; dense for bulk, scipy sparse CSR for single-cell. BAF is never stored — derive it from `Ballele` / `Tallele`.

### Final bins

Input for HATCHet3 / CalicoST. `min_snp_reads` may be a list: one job preprocesses once and writes one self-contained `MSR{msr}/` subdirectory per value.

| Mode | Location |
|---|---|
| `bulk_genotyping` | `bb_dir/MSR{msr}/{stream}/` (`stream` = `bulkWES` if the run has any `bulkWES` record, else `bulkWGS`) |
| `single_cell_genotyping` | `bb_dir/MSR{msr}/{assay_type}/` |
| `copytyping_preprocess` | `bb_dir/{assay_type}/` (not MSR-driven) |

**`bulk_genotyping`** — all bulk assays on one shared bin grid, one pseudobulk column per replicate.

| File | Contents |
|---|---|
| `bb.tsv.gz` | Bin annotations; the one grid shared by every bulk assay. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts; columns concatenate all bulk samples (per assay, normal first). |
| `bb.{depth,rdr}.npz` | Depth, and RDR for the tumor columns only. Denominator per `rdr_normalization`. |
| `sample_ids.tsv` | One row per sample, in matrix-column order. |

**`single_cell_genotyping`** — all non-bulk assays on one shared grid (one pseudobulk column per replicate x assay); the grid and sample sheet are duplicated into each per-assay subdir.

| File | Contents |
|---|---|
| `bb.tsv.gz` | The shared bin grid; identical copy in each subdir. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Per-assay allele counts (bins x cells). |
| `bb.Xcount.npz` | Per-assay native counts, same shape and column order as the allele matrices. **scATAC**: deduped fragment counts (each fragment counted once by its midpoint, from the raw fragments). **scRNA / VISIUM**: UMI counts from the h5ad, each gene assigned to its largest-overlap bin. |
| `multi_snp.tsv.gz`, `multi_snp.{T,A,B}allele.npz` | Multi-SNP diagnostic groups; MSR-independent, identical across subdirs. |
| `barcodes.tsv.gz`, `barcodes.full.tsv.gz` | Per-assay cell barcodes. |
| `sample_ids.tsv` | One row per replicate x assay; identical copy in each subdir. |

**`copytyping_preprocess`** — per assay, flat.

| File | Contents |
|---|---|
| `cnv_segments.tsv` | BB block annotations. |
| `bb.{Xcount,Tallele,Aallele,Ballele}.npz` | Per-block count matrices. |
| `barcodes.tsv.gz`, `barcodes.full.tsv.gz`, `sample_ids.tsv` | As above. |

### Intermediates

| Directory | Contents |
|---|---|
| `snp_dir/` | `chr{chrname}.vcf.gz` (bi-allelic SNPs); `pseudobulk_{modality}/cellSNP.*` and `pseudobulk_snp_statistics.tsv` (single-cell). |
| `phase_dir/` | `chr{chrname}.vcf.gz` (phased); `phased_het_snps.vcf.gz(.tbi)`; `germline_snp_statistics.tsv`; `genetic_map.tsv.gz` (eagle/shapeit). |
| `pileup_dir/` | One cellsnp-lite dir per `{assay_type}_{dataset_id}`; bulk also `{assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz`, `window.dp.npz`, `window.tsv.gz`, `depth_statistics.tsv`. |
| `allele_dir/` | Per `{stream}` (bulk) or `{assay_type}` (single-cell): `snps.tsv.gz`, `snp.{T,A,B}allele.npz`, `sample_ids.tsv`, and for single-cell `barcodes{,.full}.tsv.gz`, `unique_snp_ids.npy`. |
| `bb_dir/{assay_type}.h5ad` | Gene x cell AnnData (single-cell RNA / spatial); MSR-independent, so it sits flat. |

### TSV columns

| File | Columns |
|---|---|
| `snps.tsv.gz` | `#CHR`, `POS`, `POS0`, `START`, `END`, `GT`, `PHASE` (0 = B-allele is ALT, 1 = B-allele is REF), `region_id`, `feature_id` (`;`-joined overlapping GTF genes, `intergenic` if none), `feature_type` (exon/intron/intergenic). Bulk adds `PS` (phase set) when longphase emits it. |
| `bb.tsv.gz` | `#CHR`, `START`, `END`, `#SNPS`, `region_id`, `switchprobs`, `feature_id` (deduped union of the bin's SNP genes). `cnv_segments.tsv` carries the same `feature_id`. |
| `multi_snp.tsv.gz` | `#CHR`, `START`, `END`, `START0`, `END0`, `region_id`, `feature_id`, `feature_type`, `#SNPS`, `BLOCKSIZE`, `multi_id`, `switchprobs`. |
| `sample_ids.tsv` | `SAMPLE` (`{patient_id}_{dataset_id}`), `SAMPLE_NAME` (patient ID), `REP_ID`, `sample_type`. Under `bb_dir/` it also carries `assay_type`; row order always matches the matrix columns. |
| `germline_snp_statistics.tsv` | Per-chromosome counts: het_phased, het_unphased, hom_alt, hom_ref. |
| `depth_statistics.tsv` | Per-assay read-depth summary from `rd_correct`. |

### Barcodes (single-cell)

| File | Contents |
|---|---|
| `barcodes.tsv.gz` | One `{BARCODE}_{REP_ID}` per row, no header. |
| `barcodes.full.tsv.gz` | Two columns with header, `REP_ID` and `BARCODE`, in matrix-column order; `BARCODE` keeps the same `{BARCODE}_{REP_ID}` form. |
| `unique_snp_ids.npy` | SNP ids as `{chr}_{pos}`. |

### QC (`qc_dir/`)

One multi-page PDF per rule, flat:

| File | Contents |
|---|---|
| `phase_and_concat.{stream_or_assay}.pdf` | SNP allele frequency (unphased REF/TOTAL, then phased BAF) and SNP depth histogram. |
| `rd_correction.{assay_type}.pdf` | Read-depth bias correction (bulk): RD before/after, GC/MAP/RT diagnostics. |
| `combine_counts.{stream_or_assay}.MSR{msr}.pdf` | Binning QC, one PDF per `min_snp_reads` value. |
| `cnv_segmentation.{assay_type}.pdf` | SNP- and BB-level BAF (`copytyping_preprocess`). |
