# Universal Genotyping Pipeline (pre-release)

Universal Genotyping Pipeline inputs BAMs from various assay types (bulk LR/SR WGS/WES, scRNA-seq, scATAC-seq, VISIUM, VISIUM HD 3prime), performs panel-based SNP genotyping, population-based/read-based haplotype phasing, read counting, sequencing bias correction, and segmentation, and outputs genomic segment-level data matrices for downstream applications.

---

## Installation

Requires [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/) >= 9. Every tool dependency is declared per-rule under `workflow/envs/` (`base.yaml` — Python scientific stack; `tools.yaml` — bcftools, cellsnp-lite, mosdepth, samtools, tabix; `phase.yaml` — eagle2, shapeit5, longphase) and built by Snakemake itself:

```sh
snakemake --profile /path/to/workflow/profile/ \
    --conda-create-envs-only --cores 1 \
    -s /path/to/workflow/Snakefile
```

> **Note:** the default `conda-prefix` in `profile/config.yaml` is the relative path `.snakemake/conda`. After creating the environments, change it to an absolute path so they are reused when running with `--directory`.

---

## Quick start

Prepare a sample file with [schema](docs/sample_sheet.md) and [template](resources/templates/samples.json), then validate:

```sh
python resources/scripts/validate_sample_file.py /path/to/samples.json --check-files
```

Prepare a configuration file with [template](resources/templates/config.yaml), then run:

```sh
snakemake --profile /path/to/workflow/profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.json sample_id=<PATIENT_ID>
```

Preview with `--dry-run` (`-n`); resume a crashed run with `--rerun-incomplete`.

---

## Documentation

| Document | Description |
|----------|-------------|
| [docs/workflow.md](docs/workflow.md) | The three modes: stages, rules, and where each writes. |
| [docs/sample_sheet.md](docs/sample_sheet.md) | Sample-file schema, remote inputs, validation. |
| [docs/reference.md](docs/reference.md) | Every config key and default; every output file and column. |
| [docs/tutorials/bulk_genotyping.hg38_wgs.md](docs/tutorials/bulk_genotyping.hg38_wgs.md) | Tutorial: bulk WGS genotyping on hg38. |
| [config/config.yaml](config/config.yaml) | Default configuration (auto-loaded by the Snakefile). |
| [resources/templates/](resources/templates/) | User config and sample-file templates. |
| [resources/README.md](resources/README.md) | External resources: panels, genetic maps, window BEDs. |
