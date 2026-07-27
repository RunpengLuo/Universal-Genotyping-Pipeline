# Bulk Genotyping

This documentation covers input preparation and result intepretation for bulk genotyping using **short-read** WGS/WES and/or **long-read** (e.g., PacBio HiFi, Oxford Nanopore) sequencing data to run [HATCHet](https://github.com/raphael-group/hatchet). Refer to [Installation](installation.md) and [run.md](./run.md) for Snakemake pipeline installation and execution instructions.

<p align="center">
  <img src="imgs/rulegraph.bulk_genotyping.png" alt="bulk_genotyping rule graph" width="500">
</p>

## Table of Contents
1. [Input](#input) <br>
2. [Output](#output) <br>

## Input

### Sample file
A sample sheet in JSON format is required to specify the locations and data configurations for input datasets. Detailed JSON format can be found at [sample_sheet.md](./sample_sheet.md), Here is an example for one tumor WGS dataset `T1` with matched normal WGS sample `N1` from patient `HT001`. 

```json
{
  "version": 1,
  "samples": [
    {
      "sample_id": "HT001",
      "dataset_id": "N1",
      "assay_type": "bulkWGS",
      "sample_type": "normal",
      "files": {
        "alignment": "/data/HT001/normal.bam",
        "alignment_index": "/data/HT001/normal.bam.bai"
      }
    },
    {
      "sample_id": "HT001",
      "dataset_id": "T1",
      "rdr_base_dataset_id": "N1",
      "assay_type": "bulkWGS",
      "sample_type": "tumor",
      "files": {
        "alignment": "/data/HT001/tumor.cram",
        "alignment_index": "/data/HT001/tumor.cram.crai"
      }
    }
  ]
}
```

> [!IMPORTANT]
> Here are a few important constraints for sample files:
> - All alignment files must come from same reference version.
> - Each tuple (`sample_id`, `dataset_id`) defines a unique dataset.
> - BAM file (`alignment`) must be sorted, and its index file (`alignment_index`) must present!

> [!NOTE]
> (Experimental!) For long-read data, a pre-defined confident SV breakpoints in BEDPE format can also be provided via `files.breakpoint_bedpe` such that segmentation will avoid to over-segment across the breakpoints. See [sample_sheet.md](sample_sheet.md#files).

### Config file

A Snakemake config file is required to specify the runtime configurations. Copy the [template](../resources/templates/config.yaml) and adjust following parameters. detailed description can be found at [reference.md](reference.md#configuration).

1. specify workflow mode, assay types (`assay_types`), and patient informations.

```yaml
workflow_mode: "bulk_genotyping"
assay_types: ["bulkWGS"]

sample_id: HT001
sample_file: /path/to/samples.json
chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
```

2. specify the paths to reference files. Commonly used reference versions for human (`hg19`, `hg38`, `T2T-CHM13v2.0`) and mouse (`mm10`) have pre-built files at [../resources/data/](../resources/data/). Here is an example configuration for `hg38`. See [../resources/README.md](../resources/README.md) for detailed descriptions and public URLs.

```yaml
reference_version: hg38
reference: /path/to/reference.fasta
genome_size: resources/data/hg38.chrom.sizes
region_bed: resources/data/hg38.regions.bed
window_bed: resources/data/windows.1kbp.hg38.bed.gz
mappability_bed: resources/data/hg38.mappability.bed.gz
blacklist_bed: resources/data/hg38-blacklist.v2.bed.gz
gene_blacklist_file: resources/data/ig_gene_list.txt
gtf_file: /path/to/gencode.v38.annotation.gtf.gz
```

> [!NOTE]
> (Experimental!) When SV breakpoints are provided, `window_bed` will be ignored and re-built according to `window_size` (default is 1kbp) and breakpoints.

3. specify the targeted positions (`snp_targets`) and list of normal datasets (`genotype_dataset_ids`, default is all normal samples) for germline SNPs genotyping via bcftools. See [snp-panels](../resources/README.md#snp-panels) for details.

```yaml
snp_targets: /path/to/target_positions
```

4. our pipeline supports various haplotype phasing softwares (`phaser`) including Eagle2, Shapeit5, and LongPhase.
- For short-read phasing via Eagle2 and Shapeit5, genetic map file (`gmap_path`, see [genetic-maps](../resources/README.md#genetic-maps)) and population haplotype panel (`phasing_panel`, see [phasing-panels](../resources/README.md#phasing-panels)) are required. 
- For long-read phasing via LongPhase, genetic map and haplotype panel are ignored. Set `params_longphase.extra_params` according to specific long-read sequencing technology (e.g., `"--pb"` for Pacbio HiFi).

An example for Eagle2 as follows.
```yaml
phaser: "eagle"
phasing_panel: /path/to/1kGP_3202_hg38/phasing_panel
gmap_path: /path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
```

> [!NOTE]
> If germline (phased) Het SNPs information is already exist, user may also specify the path via `het_snp_vcf` and set `het_snp_vcf_phased` to indicate if the VCF file is phased or not. This will skip the germline SNP genotyping (and haplotype phasing if `het_snp_vcf_phased=true`).

5. By default, our pipeline performs dataset-specific read-depth bias correction (`gc_correct_method`, default is median regression, lowess is also provided) against covariates including GC-content and replication timing (RT). For long-read sequencing datasets where GC bias are less common, user can disable them by setting `gc_correct` (and `rt_correct`) to `false` as follows.

```yaml
params_count_reads:
  gc_correct: false
  rt_correct: false
```

6. The final step of the pipeline is to perform variable-length segmentation over genomic windows jointly across all tumor samples and obtain genomic bin by sample read-depth ratio (RDR), phased B-allele counts, and total-allele count matrices. each value in minimum-SNP-covering reads parameter (`min_snp_reads`) gives one segmentation result. We recommend setting `min_snp_reads` to a list of values and manually pick the lowest parameter that gives reliable BAF signals. In practice, we recommend using `MSR5000` as the final output.
```yaml
params_combine_counts:
  min_snp_reads: [500, 1000, 3000, 5000, 10000]
```

## Output

Here we show the key results and visualizations from bulk genotyping.

```text
<out_dir>/
  bb/
    MSR{msr}/                          # one subdir per min_snp_reads value
      bulk/                            # one joint grid over all bulk assays (WGS/WGS-lr/WES)
        bb.tsv.gz                      # bin annotations (grid shared by every matrix below)
        bb.{Tallele,Aallele,Ballele}.npz   # phased allele counts, bins x samples
        bb.{depth,rdr}.npz             # read depth (all samples) and RDR (tumor columns)
        sample_ids.tsv                 # one row per sample, in matrix-column order
  qc/
    rd_correction.{assay_type}.pdf     # read-depth bias correction, one per assay
    combine_counts.bulk.MSR{msr}.pdf   # binning QC, one per min_snp_reads value
```

All final bins live under `bb_dir/MSR{msr}/bulk/`, one `MSR{msr}/` subdirectory per `min_snp_reads` value. Every `.npz` is a dense matrix whose rows are the bins of `bb.tsv.gz` (same order) and whose columns are the samples of `sample_ids.tsv` (same order); BAF is never stored, derive it as `Ballele / Tallele`. Refer to [Final bins](reference.md#final-bins) for the full per-file contract.

### `bb.tsv.gz`

Bin annotations, one row per bin; the one grid shared by every `bb.*.npz`, its row order defining the matrix rows. Refer to [TSV columns](reference.md#tsv-columns) for the column definitions.

### `bb.Tallele.npz`, `bb.Aallele.npz`, `bb.Ballele.npz`

Phased allele-count matrices (bins x samples), columns concatenating all bulk samples. Inspect `qc/combine_counts.bulk.MSR{msr}.pdf` (one per `min_snp_reads`) to compare BAF signal across bin sizes.

### `bb.depth.npz`

Read depth per bin (bins x samples), all samples, aggregated from the bias-corrected window depth. Inspect `qc/rd_correction.{assay_type}.pdf` (RD before/after correction, GC/MAP/RT diagnostics) to confirm the depth is well corrected, especially for `bulkWES` whose capture-enrichment structure leaves a noisier post-correction signal.

### `bb.rdr.npz`

Read-depth ratio (RDR) for the tumor columns only (bins x tumor-samples). The denominator follows `rdr_normalization`: a tumor's `rdr_base_dataset_id` (a matched same-platform normal) when set, else the genome-wide median. Same bin rows as `bb.tsv.gz`; feeds HATCHet3 alongside `Ballele`/`Tallele`.

### `sample_ids.tsv`

One row per sample, its row order matching the columns of every `bb.*.npz`. Refer to [TSV columns](reference.md#tsv-columns) for the column definitions.
