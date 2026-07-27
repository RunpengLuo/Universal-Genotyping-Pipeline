# Pipeline Guide

This pipeline guide covers the installation, configuration, and execution of our pipeline. For specific workflow mode, see the tutorials ([bulk-genotyping](bulk_genotyping.md), [single-cell-genotyping](single_cell_genotyping.md), [copytyping-preprocess](copytyping_preprocess.md)).

## Installation

`Universal-Genotyping-Pipeline` requires a 64-bit Linux system and [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/) (version 9 or newer). Create the runner environment from [`environment.yaml`](../environment.yaml) and activate it:

```sh
conda env create -f environment.yaml
conda activate genotyping-env
```

### Pipeline Dependencies
Pipeline dependencies live under `workflow/envs/`, one file per tool group, so `--use-conda` builds only what a run needs:

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

Build the conda environments for pipeline dependencies once:

```sh
snakemake --profile profile/ \
    --conda-create-envs-only --cores 1 \
    -s workflow/Snakefile
```

## Configure

You need three things before a run: a **profile**, a **sample file**, and a **config file**.

### 1. Profile

The profile at [`profile/config.yaml`](../profile/config.yaml) holds run-wide settings. The main keys to check:

| Key | What it does |
|-----|--------------|
| `cores` | Max CPU cores the run may use. |
| `use-conda` | Keep `true` so per-rule tool environments are used. |
| `conda-prefix` | Where the tool environments live; set to an **absolute path** (see the note above). |
| `resources: downloads` | How many remote (URL) input files download at once; default `2`. |
| `local-storage-prefix` | Where remote files are downloaded; point it at a large scratch disk if inputs are URLs. |

For local BAM inputs the defaults are fine.

> [!IMPORTANT]
> The default `conda-prefix` is the relative path `.snakemake/conda` w.r.t. current working directory. After the environments are built, change it to an absolute path to avoid re-building the environment for every new run.

### 2. Sample file

The sample file lists your input datasets. Copy the [template](../resources/templates/samples.json), fill it in using the [schema](sample_sheet.md), then check it:

```sh
python resources/scripts/validate_sample_file.py /path/to/samples.json --check-files
```

`--check-files` confirms every local path exists. See the per-mode tutorials for full examples: [bulk](bulk_genotyping.md), [single-cell](single_cell_genotyping.md), [copytyping](copytyping_preprocess.md).

### 3. Config file

The config file sets the workflow mode, reference files, and parameters. Copy the [template](../resources/templates/config.yaml) and modify it; every key is described in [reference.md](reference.md#configuration).

## Running the pipeline

Run with your profile, Snakefile, config file, output directory, and sample file:

```sh
snakemake --profile profile/ \
    -s workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.json sample_id=<PATIENT_ID>
```

`--config sample_file=... sample_id=...` picks which sample file and which patient to run; they can also be set inside the config file instead.

> [!TIP]
> Always preview first with a **dry run** - add `--dry-run` (`-n`). It lists the jobs Snakemake would run without running any, so you can confirm the plan and catch config or sample-file mistakes early:
> ```sh
> snakemake --profile profile/ -s workflow/Snakefile \
>     --configfile /path/to/my_config.yaml --directory <output_dir> \
>     --config sample_file=/path/to/samples.json sample_id=<PATIENT_ID> \
>     -n
> ```

If a run stops partway, resume it with `--rerun-incomplete`. Outputs land under `<output_dir>`; see each mode's tutorial for what to expect.
