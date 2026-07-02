# Tutorial: Bulk WGS Genotyping (hg38)

This tutorial walks through running the `bulk_genotyping` pipeline on paired normal/tumor bulk WGS samples aligned to GRCh38.

## 1. Prerequisites

Install [conda](https://github.com/conda-forge/miniforge) and [Snakemake](https://snakemake.readthedocs.io/) (>= 9). We recommend using the [libmamba solver](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community) for faster dependency resolution (`conda config --set solver libmamba`).

## 2. Sample sheet

Create a TSV file (e.g., `samples.tsv`) with one row per BAM. For bulk WGS you need five columns:

```tsv
SAMPLE	REP_ID	assay_type	sample_type	PATH_to_bam
HT001	N1	bulkWGS	normal	/data/HT001/normal.bam
HT001	T1	bulkWGS	tumor	/data/HT001/tumor.bam
```

You may specify data from multiple patients in same sample sheet, but the workflow will only process one patient at a time depends on specified `sample_id` in config file. See [docs/reference.md](../reference.md) for detailed format.

## 3. Config

Copy `resources/templates/config.yaml` and edit the sections below. Defaults are auto-loaded from `config/config.yaml`; only the fields that need changing are shown.

### Workflow mode

```yaml
workflow_mode: "bulk_genotyping"
assay_types: ["bulkWGS"]
```

### Input

```yaml
sample_id: HT001
sample_file: /path/to/samples.tsv
chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
```

### Reference

| Resource | Description |
|----------|-------------|
| `windows.1kbp.hg38.bed.gz` | Pre-built 1 kb window BED with GC, mappability, and replication timing (included in `resources/data/`) |
| `hg38-blacklist.v2.bed.gz` | ENCODE blacklist v2 (included in `resources/data/`) |
| `hg38.regions.bed` | Chromosome region definitions, excluded centromeric region and acrocentric short arms (13/14/15/21/22 p-arm) (included in `resources/data/`) |
| **GTF annotation** | GENCODE v38 ([download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)) |

```yaml
reference_version: hg38
reference: /path/to/reference.fasta
genome_size: resources/data/hg38.regions.bed
region_bed: resources/data/hg38.regions.bed
window_bed: resources/data/windows.1kbp.hg38.bed.gz
blacklist_bed: resources/data/hg38-blacklist.v2.bed.gz
gtf_file: /path/to/gencode.v38.annotation.gtf.gz
```

### SNP and Phasing panels
We recommend the [1kGP n=3,202 high-coverage](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/) panel. Download and prepare the 1kGP n=3,202 panel (requires bcftools, bgzip, tabix):

```bash
bash resources/scripts/process_1kGP_3202_panel.sh --ref hg38 /path/to/1kGP_3202
```

This script will produce three directories under the output path:

```yaml
snp_panel: /path/to/snps.vcf.gz
snp_targets: /path/to/target_positions
phasing_panel: /path/to/phasing_panel
```

### Phasing
| Resource | Description |
|----------|-------------|
| **Eagle2** | Phasing tool ([download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz)) — includes genetic maps in `tables/` |
| **SHAPEIT5** | Phasing tool ([GitHub](https://github.com/odelaneau/shapeit5)) — includes genetic maps in `resources/maps/` |

Pick one of Eagle2 or SHAPEIT5. Both require a genetic map (gmap) bundled with the tool distribution. SHAPEIT5 is recommended.

One config key, `gmap_path`, locates the gmap. Use `{chrname}` as a placeholder for SHAPEIT5 (per-chrom files); for Eagle2 (single file) use a literal path.

**Eagle2** — gmap shipped at `{Eagle2_dir}/tables/genetic_map_hg38_withX.txt.gz`.

```yaml
phaser: "eagle"
gmap_path: "/path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
```

**SHAPEIT5** — gmap shipped at `{shapeit5_dir}/resources/maps/b38/chr{N}.b38.gmap.gz`.

```yaml
phaser: "shapeit"
gmap_path: "/path/to/shapeit5/resources/maps/b38/chr{chrname}.b38.gmap.gz"
```

### Count reads
By default, we apply quadratic median regression to correct GC content and replication timing biases while being robust to copy-number signals. Use read-depth plots and KDE scatter plot for diagnostic purpose under `<qc_dir>/rd_correction.rd_correct.<assay_type>.*.pdf`.

### Combine counts

`min_snp_reads` and `min_snp_per_block` balance the bin size against phase-switch errors from reference phasing. `min_snp_reads` accepts a **list**: one `combine_counts` job then loads/preprocesses once and writes one bin grid per value under `bb_dir/MSR{msr}/{stream}/`, with a QC PDF per value at `<qc_dir>/combine_counts.{stream}.MSR{msr}.pdf`. Compare the bb BAF/RDR genome plots across `MSR{msr}` and keep the subdir whose bins are neither too sparse nor too noisy for your downstream tool.

```yaml
params_combine_counts:
  min_snp_reads: [500, 1000]
  min_snp_per_block: 10
```

## 4. Profile and running
You may use the profile at `profile/config.yaml` to adjust snakemake parameters like `cores` to match your machine.

### Pre-create conda environments
```bash
snakemake --profile /path/to/profile/ \
    --conda-create-envs-only --cores 1 \
    -s /path/to/workflow/Snakefile
```

> **Note:** The default `conda-prefix` in `profile/config.yaml` is the relative path `.snakemake/conda`. After creating the environments, change it to an absolute path (e.g., `/path/to/workflow/.snakemake/conda`) so that conda environments are reused correctly when running with `--directory`.

### Run the pipeline

```bash
snakemake --profile /path/to/profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory /path/to/output_dir \
    --config sample_file=/path/to/samples.tsv sample_id=HT001
```

Defaults are auto-loaded from `config/config.yaml`. The `--configfile` and `--config` flags override specific values, so you can reuse the same config for different patients.

## 5. Results

Final outputs required by HATCHet3 can be found in `<out_dir>/<bb_dir>/MSR{msr}/{stream}/` (`stream` = `bulkWGS` for the WGS family or `bulkWES`; e.g. `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele,depth,rdr}.npz`). See [docs/reference.md](../reference.md) for the full output specification including intermediate files.

## 6. Rerunning from intermediate results

Snakemake automatically skips rules whose outputs already exist. To resume after a crash, add `--rerun-incomplete` to rerun rules that left partial outputs. Adding `--keep-incomplete` to the initial run preserves partial outputs on failure instead of deleting them.

```bash
snakemake --profile /path/to/profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory /path/to/output_dir \
    --config sample_file=/path/to/samples.tsv sample_id=HT001 \
    --rerun-incomplete
```

Add `--dry-run` (`-n`) to preview what will be rerun before executing.
