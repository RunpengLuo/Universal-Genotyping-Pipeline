# Bulk Genotyping

This documentation covers input preparation and result intepretation for preprocessing **short-read** WGS/WES and/or **long-read** (e.g., PacBio HiFi, Oxford Nanopore) sequencing data to run [HATCHet](https://github.com/raphael-group/hatchet). Refer to the [README](../README.md) for Snakemake pipeline installation and execution instructions.

## Table of Contents
1. [Overview](#overview) <br>
2. [Input](#input) <br>
3. [Output](#output) <br>

## Overview

The rule graph below shows the stages of the bulk genotyping workflow.

<p align="center">
  <img src="imgs/rulegraph.bulk_genotyping.png" alt="bulk_genotyping rule graph" width="500">
</p>

### Features
- Joint genotyping and segmentation across multiple bulk samples of mixed **short-reads** and **long-reads** DNA sequencing platforms on one shared genomic bin coordinates.
- Large data hosted on cloud or FTP servers can be streamed rather than download to local storage. (`remote_mode=stream`).
- Sequencing read-depth bias correction (GC-content, mappability, replication timing).
- Multi-sample adaptive genomic bin segmentation to achieve sufficient count statistics.

### Limitations
- We require a matched-normal sample for genotyping germline SNPs. We plan to extend the pipeline to genotype germline SNPs using tumor-only samples in the future.
- Automated model-selection for segmentation parameters with respect to sequencing coverages and segmentation variances.

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
> - All alignment files must come from same reference version.
> - Each tuple (`sample_id`, `dataset_id`) defines a unique dataset.
> - BAM file (`alignment`) must be sorted, and its index file (`alignment_index`) must present!

> [!NOTE]
> (Experimental!) For long-read data, a pre-defined confident SV breakpoints in [BEDPE](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format) format can also be provided via `files.breakpoint_bedpe` such that segmentation will avoid to over-segment across the breakpoints. See [sample_sheet.md](sample_sheet.md#files).

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
species: human
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
> (Experimental!) When SV breakpoints are provided, `window_bed` will be ignored and re-built according to `window_size` (default: 1kbp) and breakpoints.

3. specify the targeted positions (`snp_targets`) and list of normal datasets (`genotype_dataset_ids`, default is all normal samples) for germline SNPs genotyping via [bcftools](https://github.com/samtools/bcftools). See [snp-panels](../resources/README.md#snp-panels) for details.

```yaml
snp_targets: /path/to/target_positions
```

> [!TIP]
> - If a set of confident germline (phased) Het SNPs information already exist, user may specify the path via `het_snp_vcf` and set `het_snp_vcf_phased` to indicate if the VCF file is phased or not. This will skip the germline SNP genotyping (and haplotype phasing if `het_snp_vcf_phased=true`).
> - For long-read datasets, set `params_bcftools.extra_params` to the matching bcftools mpileup platform preset so genotyping and pileup use the correct long-read error model: `-X ont-sup` (Oxford Nanopore) or `-X pacbio-ccs` (PacBio HiFi). Run `bcftools mpileup -X list` for all available profiles.

4. our pipeline supports various haplotype phasing softwares:
- For short-read phasing via [Eagle2](https://github.com/poruloh/Eagle) (preferred) and [Shapeit5](https://github.com/odelaneau/shapeit), genetic map file (`gmap_path`, see [genetic-maps](../resources/README.md#genetic-maps)) and population haplotype panel (`phasing_panel`, see [population-haplotype-panels](../resources/README.md#population-haplotype-panels)) are required. 
- For long-read phasing via [LongPhase](https://github.com/twolinin/longphase), genetic map and haplotype panel are ignored, all matched-normal long-read data will be used as input. Set `params_longphase.extra_params` according to specific long-read sequencing technology (e.g., `"--pb"` for Pacbio HiFi).

Here is an example for Eagle2:
```yaml
phaser: "eagle"
phasing_panel: /path/to/1kGP_3202_hg38/phasing_panel
gmap_path: /path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz
```

> [!TIP]
> For long-read phasing via LongPhase, set `params_longphase.extra_params` to the sequencing platform flag: `--ont` (Oxford Nanopore) or `--pb` (PacBio HiFi/CCS).

5. For each dataset, our pipeline (by default) performs read-depth sequencing bias correction (`gc_correct_method`, default: median regression) against covariates including GC-content and replication timing (RT). User may disable them:

```yaml
params_count_reads:
  gc_correct: false
  rt_correct: false
```

> [!NOTE]
> We recommend user to apply this correction and inspect the effects via `<qc_dir>/rd_correction.{assay_type}.pdf`.

6. The final step of the pipeline is to perform variable-length segmentation over genomic windows jointly across all tumor samples and obtain genomic bin by sample read-depth ratio (RDR), phased B-allele counts, and total-allele count matrices. each value in minimum-SNP-covering reads parameter (`min_snp_reads`) gives one segmentation result. We recommend setting `min_snp_reads` to a list of values and inspect the QC plots for varying `min_snp_reads`, then pick the lowest value that gives reliable BAF signals.
```yaml
params_combine_counts:
  min_snp_reads: [100, 500, 1000, 2000, 3000, 5000, 7500, 10000]
  max_blocksize: 500000 # default: 0.5MB.
```

> [!TIP]
> 1. For high-coverage (>=30x) data, we recommend to use `min_snp_reads>1000`. In practice we usually use 5000.
> 2. For low-coverage/targeted data, we recommend to use `min_snp_reads<200`. In practice we usually use 100.

## Output

Refer to [Final bins](reference.md#final-bins) for the full specification of each file:

```text
<out_dir>/
  ...
  bb/
    MSR{msr}/                          # one subdir per min_snp_reads value
      bulk/                            # one joint grid over all bulk assays (WGS/WGS-lr/WES)
        bb.tsv.gz                      # bin annotations (grid shared by every matrix below)
        bb.{Tallele,Aallele,Ballele}.npz   # phased allele counts, bins x samples
        bb.{depth,rdr}.npz             # read depth (all samples) and RDR (tumor columns)
        sample_ids.tsv                 # one row per sample, in matrix-column order
  qc/
    genotype_snp_qc.pdf                # het vs hom-alt ref-AF diagnostic (qc_genotype_snps)
    phase_and_concat.bulk.pdf          # SNP allele frequency + per-dataset depth (phase_and_concat)
    rd_correction.{assay_type}.pdf     # read-depth bias correction, one per assay
    combine_counts.bulk.MSR{msr}.pdf   # binning QC (segmentation, genome-wide RDR/BAF, RDR-vs-BAF 2D), one per min_snp_reads value
```
