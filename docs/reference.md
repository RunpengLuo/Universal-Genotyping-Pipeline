# Reference

## Table of Contents
- [Dependencies](#dependencies)
- [Sample File](#sample-file)
- [Configuration](#configuration)
  - [Input Data](#input-data)
  - [Parameters](#parameters)
- [Outputs](#outputs)
  - [Genomic unit levels](#genomic-unit-levels)
  - [Final bins](#final-bins)
  - [Intermediates](#intermediates)
  - [TSV columns](#tsv-columns)
  - [Barcodes (single-cell)](#barcodes-single-cell)
  - [QC](#qc-qc_dir)

---

## Dependencies

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

## Sample File
Refer to spec **[sample_sheet.md](sample_sheet.md)** and [templates](../resources/templates/).

---

## Configuration
Defaults in `config/config.yaml`, template in [templates](../resources/templates/). Override with `--config key=value`.

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
| `snp_panel` | Genotyping | Population SNP VCF. |
| `snp_targets` | Bulk genotyping | Per-chromosome position files; build via `resources/scripts/build_snp_targets.sh`. |
| `phaser` | Genotyping | `eagle` \| `shapeit` \| `longphase`. |
| `phasing_panel` | eagle/shapeit | Per-chromosome BCF reference panel directory. |
| `gmap_path` | eagle/shapeit | Genetic map; `{chrname}` placeholder for per-chromosome maps (SHAPEIT5), literal path for a single map (Eagle2). |
| `genotype_dataset_ids` | Optional | `dataset_id`s piled up to call germline SNPs. Empty -> auto (normal before tumor, short-read before long-read). >1 are pooled in one `mpileup` and must share an `@RG SM` tag. |
| `phase_dataset_ids` | Optional | datasets used for long-read phasing inputs |
| `het_snp_vcf` | Optional; required for `copytyping_preprocess` | Pre-computed gHET VCF. |
| `het_snp_vcf_phased` | Optional | Default `true`: the input `het_snp_vcf` is phased or not. |
| `bb_file` | copytyping_preprocess | Pre-computed bb annotations TSV. |

> [!NOTE] `genome_size` defines the chromosome names of `reference` and input alignments files. Internally, the workflow detects the chr-notation of input data and normalize the final output files with chr-prefix notation.

### Parameters

Defaults are those in `config/config.yaml`.

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

| Field | Default | Description |
|---|---|---|
| `window_size` | `1000` | Window size (bp); fixed tiling of the segment BED, shared by every assay of the run. |

#### `params_bcftools`
Used by `genotype_snps_bulk` and `pileup_snps_bulk_bcftools` (bulk het-SNP read counting).

| Field | Default | Description |
|---|---|---|
| `min_mapq` | `20` | Skip alignments below this mapping quality (genotype + pileup). |
| `min_baseq` | `20` | Skip bases below this base quality (genotype + pileup). |
| `min_dp` | `5` | Minimum depth to keep a site (genotype only). |
| `max_depth` | `1000` | Per-file depth cap in `mpileup` (genotype + pileup). |
| `min_qual` | `30` | Minimum variant QUAL (genotype only). |
| `extra_params` | `""` | Extra `mpileup` flags applied to bulk genotyping + pileup; use a read-type preset, e.g. `-X ont-sup` (ONT), `-X pacbio-ccs`, or `""`/`-X illumina` for short-read. |

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
| `rt_correct` | `true` | Also correct replication-timing bias. Gates the ENCODE Repli-seq fetch: `false` skips the 16 bigWigs and the liftOver, and the window BED gets no `REPLI` column. |
| `samplesize` | `50000` | Max ideal bins used to fit the correction curve. |
| `routlier` | `0.01` | Upper quantile of read count dropped as outlier. |
| `doutlier` | `0.001` | Top/bottom quantile of the GC/mappability domain dropped as outlier. |
| `min_mappability` | `0.9` | Bins below this mappability are excluded from the correction fit AND set to NaN in `window.dp.npz`, for either `gc_correct_method` and even when `gc_correct: false`. |

> [!NOTE]
> One `rd_correct` job corrects every bulk dataset onto the one window grid, so
> `pileup/bulk/window.dp.npz` is `(windows x bulk datasets)` with no row dropped. A bin with
> zero coverage corrects to `0.0` - zero is an observation, and a deleted region keeps its
> signal - and NaN marks only a correction the model cannot define. With the default
> `gc_correct_method: median` that is exactly one case, and it applies to every dataset
> alike: mappability below `min_mappability`. `lowess` adds two dataset-specific cases, since
> it fits per dataset and does not extrapolate: a bin whose GC or mappability falls outside
> the fitted range, and (with `rt_correct: true`) a bin whose `REPLI` is missing from the
> lifted-over track. `combine_counts` masks NaN per dataset column when it aggregates windows
> into bbs.

#### `params_combine_counts`
Used by `combine_counts` (bulk) and `combine_counts_nonbulk` (single-cell).

| Field | Default | Description |
|---|---|---|
| `min_snp_reads` | `[500, 1000]` | SNP-read target per bin, applied to every tumor column in both modes (bulk: one per sample, WGS/WGS-lr/WES alike; single-cell: one pseudobulk per tumor `dataset_id x assay_type`). A list sweeps binning: one job preprocesses once and writes one `MSR{msr}/` subdir per value. |
| `min_snp_per_bin` | `1` | Minimum SNPs per bin. |
| `gene_aware_binning` | `true` | `true` = a bin grows by whole genes (GTF `feature_id`) and never splits one; `false` = bin/SNP-granular binning. |
| `nu` | `1` | Scale of the Haldane map `(1 - exp(-2*nu*d)) / 2` turning cM distance into a switch probability. |
| `min_switchprob` | `1e-6` | Floor on that switch probability. |
| `switchprob_ps` | `5e-2` | Switch probability within one phase set (`PS`); across phase sets it is ~0.5. |
| `nsnp_multi` | `2` | SNPs per multi-SNP diagnostic group. Written in both modes, outside the `MSR{msr}/` layer. |
| `max_blocksize` | `500000` | Force a bin cut past this span, in bp (bulk; `0` = no cap). |
| `rdr_normalization` | `auto` | RDR denominator of a bulk tumor: `auto` = its `rdr_base_dataset_id` when set, else median; `median` = always median; `normal` = always `rdr_base_dataset_id`, error if any tumor lacks one. Resolved at parse time; a tumor's base should be a same-platform normal so RDR cancels platform bias. |
| `phase_flip_test` | `true` | Test each phase cluster for a haplotype flip and split it (bulk). |
| `phase_flip_epsilon` | `0.05` | Effect size of that test (bulk). |
| `phase_flip_alpha` | `0.05` | Significance level of that test (bulk). |

> [!NOTE]
> High RDR is never clipped: a focal amplification is signal, and unreliable sequence is
> already excluded by `min_mappability` and `blacklist_bed`.

> [!NOTE]
> Read target: `min_snp_reads` is applied per tumor column and a bin closes only when every
> tumor column meets it. Normals are excluded in both modes - a normal is an RDR base or a
> reference and its BAF is uninformative, so a shallow normal must not drive bin size. Bulk
> handles WES identically to WGS -- same `segment.bed`/`window_size` tiling, no
> capture-target file -- so WGS/WGS-lr/WES of one patient bin together on one bb set. WES
> depth carries capture-enrichment structure, so its post-correction RDR is noisier than WGS;
> inspect its panels in `qc/rd_correction.bulk.pdf`.

> [!NOTE]
> A bulk bb is dropped when any column carries NaN in BAF, depth or RDR, so the written
> matrices are complete. The sources:
> - **BAF** is `Ballele / Tallele`, undefined where a column's total is zero: the bb holds no
>   SNP (a whole cluster without one), or one sample has no read over the bb's SNPs (a
>   shallow dataset, a WES off-target stretch, a homozygous deletion). `min_snp_reads` guards
>   the tumor columns only, so a normal reaching zero is the usual case.
> - **Depth** is the length-weighted mean over the bb's finite windows, undefined where a
>   column has none: every window below `min_mappability`, or (under `lowess`) every window
>   outside that dataset's fitted covariate range. On the multi-SNP path a group can also
>   hold no window at all, since its groups are SNP positions and windows attach by midpoint.
> - **RDR** inherits every depth NaN and adds two: a matched-normal base whose depth is zero
>   or NaN, which removes the bb for every tumor pointing at it; and a tumor with no positive
>   finite bb anywhere, whose entire column is NaN.
>
> Single-cell writes every bb. It has no depth or RDR, and a SNP-free bb still carries
> `Xcount`.

#### `threads`
Used by all multi-thread rules.

| Field | Default | Description |
|---|---|---|
| `genotype` | `4` | Threads for genotyping. |
| `phase` | `4` | Threads for phasing. |
| `pileup` | `8` | Threads for the pileup step (bulk `bcftools mpileup`; single-cell cellsnp-lite). |
| `mosdepth` | `4` | Threads for mosdepth. |

---

## Outputs

Directories (`snp_dir`, `phase_dir`, `pileup_dir`, `allele_dir`, `bb_dir`, `qc_dir`, `log_dir`, `aux_dir`, `bench_dir`) are set in `config.yaml`, relative to `snakemake --directory`. Each rule logs to `log_dir/{rule}/...` and writes a Snakemake `benchmark:` TSV (wall time, `max_rss`, `max_vms`, `cpu_time`, ...) to `bench_dir/{rule}/...`.

`.npz` are matrices: features (SNPs or bbs) x observations (samples or cells); dense for bulk, scipy sparse CSR for single-cell. BAF is never stored — derive it from `Ballele` / `Tallele`.

### Genomic unit levels

Bulk binning nests four levels, coarse to fine. The SNP is the separate allele axis assigned into bins and bbs. The vocabulary and the naming rules behind it are in [DEVELOPER.md](DEVELOPER.md).

| Level | id | Description |
|---|---|---|
| Region (chromosome arm) | `region_id` | Derived from each `region_bed` row's coordinates, `chr1:0-121700000`. Carried for RDR/QC; a bb belongs to one arm. |
| Segment | `seg_id` | Derived from each `segment_bed` row's coordinates and stamped by `build_segment_bed`; blacklist pieces of one segment keep it. The hard bb boundary: binning clusters by `seg_id` and never merges across it. With `segment_bed` unset there is one segment per arm and `seg_id == region_id`. |
| Bin (fixed tile) | `bin_id` | Fixed tile (`window_size` = 1 kb) shared by every assay of the run, and the unit `build_adaptive_bins` merges in both modes. Bulk additionally has `run_mosdepth` count it per dataset and one `rd_correct` job bias-correct every bulk dataset onto it. Intermediate, never written out; carries `#CHR/START/END/region_id/seg_id` + `GC/MAP/REPLI`. Its file is the **window BED** (`window_bed`, `aux/windows.bed.gz`) — a legacy name for the same thing. |
| bb (merged bin) | `bb_id` | Final unit: `build_adaptive_bins` merges consecutive fixed bins within one `seg_id` until every tumor observation meets `min_snp_reads` and `min_snp_per_bin`; `max_blocksize` caps only once the read target is met. Each assay's fixed-bin depth is aggregated onto these bbs. Features of `bb.tsv.gz` and every `bb.*.npz`. |
| SNP | entry of `snps.tsv.gz` | Phased het SNP; the allele axis. Assigned to its containing bin and bb; `bb.{T,A,B}allele` aggregate a bb's SNP counts. |

### Final bins

Input for HATCHet3 / CalicoST. `min_snp_reads` may be a list: one job preprocesses once and writes one self-contained `MSR{msr}/` subdirectory per value.

| Mode | Location |
|---|---|
| `bulk_genotyping` | `bb_dir/MSR{msr}/bulk/` (one joint bb set over all bulk assays: WGS/WGS-lr/WES) |
| `single_cell_genotyping` | `bb_dir/MSR{msr}/{assay_type}/` |
| `copytyping_preprocess` | `bb_dir/{assay_type}/` (not MSR-driven) |

**`bulk_genotyping`** — all bulk assays on one shared set of bbs, one pseudobulk observation per replicate.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations; the one set shared by every bulk assay. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Allele counts; columns concatenate all bulk samples, in `sample_ids.tsv` order. |
| `bb.{depth,rdr}.npz` | Depth, and RDR for the tumor columns only. Denominator per `rdr_normalization`. |
| `sample_ids.tsv` | One row per sample, in matrix-column order. |
| `multi_snp/bulk/` | The same seven files over multi-SNP groups (`nsnp_multi` SNPs each): binning-independent, so the directory sits outside the `MSR{msr}/` layer and is written once. |

**`single_cell_genotyping`** — all non-bulk assays on one shared set of bbs (one pseudobulk observation per replicate x assay); the bbs and sample sheet are duplicated into each per-assay subdir.

| File | Contents |
|---|---|
| `bb.tsv.gz` | The shared bbs; identical copy in each subdir. |
| `bb.{Tallele,Aallele,Ballele}.npz` | Per-assay allele counts (bins x cells). |
| `bb.Xcount.npz` | Per-assay native counts, same shape and column order as the allele matrices. **scATAC**: deduped fragment counts (each fragment counted once by its midpoint, from the raw fragments). **scRNA / VISIUM**: UMI counts from the h5ad, each gene assigned to its largest-overlap bin. |
| `multi_snp/{assay_type}/bb.tsv.gz`, `multi_snp/{assay_type}/bb.{T,A,B}allele.npz` | Multi-SNP diagnostic groups (`nsnp_multi` SNPs each): binning-independent, so they sit outside the `MSR{msr}/` layer and are written once per assay, with the same filenames and schema as a bb set. |
| `barcodes.tsv.gz` | This assay's cells, the matrix column axis, in column order. |
| `sample_ids.tsv` | This assay's datasets only, one row each (a roster, not column-aligned). |

**`copytyping_preprocess`** — per assay, flat.

| File | Contents |
|---|---|
| `bb.tsv.gz` | bb annotations: every column of the given `bb_file`, plus `bb_id` and (RNA) `#feature`, with `#SNPS` and `feature_id` re-derived from this assay's SNPs. Nothing is binned here, so `switchprobs` is whatever the input carried. |
| `bb.{Xcount,Tallele,Aallele,Ballele}.npz` | Per-bb count matrices. |
| `barcodes.tsv.gz`, `sample_ids.tsv` | As above. |

### Intermediates

| Directory | Contents |
|---|---|
| `snp_dir/` | `chr{chrname}.vcf.gz` (bi-allelic SNPs); `pseudobulk_{modality}/cellSNP.*` (single-cell). |
| `phase_dir/` | `chr{chrname}.vcf.gz` (phased); `phased_het_snps.vcf.gz(.tbi)`; `germline_snp_statistics.tsv`; `genetic_map.tsv.gz` (eagle/shapeit). |
| `pileup_dir/` | One pileup dir per `{assay_type}_{dataset_id}`: bulk `bcftools.counts.tsv.gz` (bcftools REF/ALT depths at the phased loci, consumed by `phase_and_concat_bulk`), single-cell `cellSNP.*` (cellsnp-lite); bulk also `{assay_type}/out_mosdepth/{dataset_id}.regions.bed.gz` (per-dataset mosdepth) and, once for the whole run, `bulk/window.dp.npz` (corrected depth, the shared windows x every bulk dataset) and `bulk/depth_statistics.tsv`. |
| `allele_dir/` | Flat, and the same names in every mode: `snps.tsv.gz` (the SNP grid = matrix rows), `snp.{T,A,B}allele.npz` (ONE matrix), `sample_ids.tsv`, and for single-cell `barcodes.tsv.gz`. The single matrix spans every assay of the run: bulk columns are samples, single-cell columns are cells over all assays. |
| `bb_dir/{assay_type}.h5ad` | Gene x cell AnnData (single-cell RNA / spatial); MSR-independent, so it sits flat. |
| `aux_dir/` | `segment.bed` (region_id arm + seg_id segment, built from `segment_bed`; every mode). Bulk window build: `windows.bed.gz` (the one shared window BED), and `repliseq/{name}.{reference_version}.bedGraph` (Repli-seq tracks, cached across window rebuilds) when `do_repliseq` (bulk, a Repli-seq build, and `rt_correct`). Source bigWigs are hg19: an hg19 run uses them directly; an hg38 or chm13v2 run lifts them over first. |

### TSV columns

| File | Columns |
|---|---|
| `snps.tsv.gz` | `#CHR`, `POS`, `POS0`, `START`, `END`, `GT`, `PHASE` (0 = B-allele is ALT, 1 = B-allele is REF), `region_id`, `seg_id` (the segment, when the segment BED carries it), `feature_id` (`;`-joined overlapping GTF genes, `intergenic` if none), `feature_type` (exon/intron/intergenic). Bulk adds `PS` (phase set) when longphase emits it. |
| `bb.tsv.gz` | `#CHR`, `START`, `END`, `#SNPS`, `region_id`, `switchprobs`, `feature_id` (deduped union of the bin's SNP genes). `copytyping_preprocess` writes its own `bb.tsv.gz`; see that mode's table. |
| `multi_snp/{bulk_or_assay}/bb.tsv.gz` | The `bb.tsv.gz` columns above, one row per multi-SNP group. |
| `sample_ids.tsv` | `SAMPLE` (derived) then the sample-file record keys verbatim: `sample_id`, `dataset_id`, `sample_type`, `assay_type`, and `rdr_base_dataset_id` (bulk only, omitted when no tumor has one). **Bulk**: one row per matrix column, in column order; `SAMPLE` is `{sample_id}_{dataset_id}`. **Single-cell**: a roster, one row per `(dataset_id, assay_type)` and NOT column-aligned - the columns are cells, listed by `barcodes.tsv.gz`. A multiome pair shares one `dataset_id`, so there `SAMPLE` is `{sample_id}_{dataset_id}_{assay_type}` to stay unique. |
| `germline_snp_statistics.tsv` | Per-chromosome counts: het_phased, het_unphased, hom_alt, hom_ref. |
| `depth_statistics.tsv` | Per-dataset read-depth summary from `rd_correct`, all bulk datasets in one table. |

### Barcodes (single-cell)

| File | Contents |
|---|---|
| `barcodes.tsv.gz` | The single-cell column axis: one `{BARCODE}_{dataset_id}_{assay_type}` per row, no header, in matrix-column order (assay-major). The parse is positional - no assay_type contains `_` and the raw barcode may not either (asserted at write) - so the last `_` ends the key and the first ends the barcode, leaving the dataset_id free to hold underscores. |

### QC (`qc_dir/`)

One multi-page PDF per rule, flat:

| File | Contents |
|---|---|
| `phase_and_concat.{bulk_or_assay}.pdf` | SNP allele frequency (unphased REF/TOTAL, then phased BAF) and SNP depth histogram. |
| `rd_correction.bulk.pdf` | Read-depth bias correction (bulk): RD before/after, GC/MAP/RT diagnostics, one panel per dataset. |
| `combine_counts.{bulk_or_assay}.MSR{msr}.pdf` | Binning QC, one PDF per `min_snp_reads` value. |
| `combine_counts_fixed_bins.{assay_type}.pdf` | SNP- and BB-level BAF (`copytyping_preprocess`). |
