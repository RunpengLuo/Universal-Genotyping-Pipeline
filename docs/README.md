# Step-by-Step Guide

## Installation

Requires [conda](https://docs.conda.io/en/latest/) and [Snakemake](https://snakemake.readthedocs.io/) >= 9.

All tool dependencies are managed via conda environments under `workflow/envs/`:

| Env file | Tools |
|----------|-------|
| `base.yaml` | Python scientific stack (scipy, numpy, pandas, numba, anndata, scanpy, etc.) |
| `tools.yaml` | bcftools, cellsnp-lite, mosdepth, samtools, tabix |
| `phase.yaml` | eagle2, shapeit5, longphase, bcftools, tabix |


```sh
snakemake --profile /path/to/workflow/profile/ \
    --conda-create-envs-only --cores 1 \
    -s /path/to/workflow/Snakefile
```

> **Note:** The default `conda-prefix` in `profile/config.yaml` is the relative path `.snakemake/conda`. After creating the environments, change it to an absolute path (e.g., `/path/to/workflow/.snakemake/conda`) so that conda environments are reused correctly when running with `--directory`.

---

## Running the pipeline

```sh
snakemake --profile /path/to/workflow/profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory <output_dir> \
    --config sample_file=/path/to/samples.tsv sample_id=<PATIENT_ID>
```

Defaults are auto-loaded from `config/config.yaml`. Pass `--configfile` with your run-specific config (see `resources/templates/config.yaml`) to override paths and settings. `--config` flags override individual keys. `sample_id` must match a `SAMPLE` value in the sample sheet.

---

## Mode 1: `bulk_genotyping`

Bulk WGS/WES. Genotypes SNPs (bcftools), phases, computes allele counts and bias-corrected RDR.

1. Sample sheet with `assay_type` = `bulkWGS`, `bulkWGS-lr`, or `bulkWES`, including normal and tumor.
2. Set `workflow_mode: bulk_genotyping` in config.
3. **WGS:** Set `window_bed` to a pre-filtered window BED (e.g., build using `resources/scripts/build_wgs_window_bed.py` or use pre-built `resources/data/windows.1kbp.hg38.bed.gz`).
4. **WES:** Set `window_bed` to a WES window BED built using `resources/scripts/build_wes_window_bed.py` (requires `--wes_targets_bed` pointing to vendor capture targets).
5. Outputs in `bb_dir/{stream}/` (`stream` = `bulkWGS` for the WGS family or `bulkWES`):
   - `bb.tsv.gz` — bin annotations (one shared grid for all bulk assays).
   - `bb.{Tallele,Aallele,Ballele,depth,rdr}.npz` — allele, depth, and RDR matrices; columns concatenate all bulk samples (per assay, normal first), `rdr` holds the tumor columns. RDR is normalized per assay against that assay's own normal. (BAF is not stored — derive from `Ballele`/`Tallele`.)
   - `sample_ids.tsv` — sample metadata, with an `assay_type` column; row order matches the matrix columns.

---

## Mode 2: `single_cell_genotyping`

scRNA, scATAC, VISIUM, or VISIUM3prime. Pseudobulk genotyping via cellsnp-lite, per-cell pileup, binned allele counts.

1. Sample sheet with non-bulk samples. Multiome: same `REP_ID` for scRNA + scATAC. Fill in `PATH_to_barcodes` and `PATH_to_10x_ranger`.
2. Set `workflow_mode: single_cell_genotyping` in config.
3. All of the sample's non-bulk assays are jointly segmented on **one shared bin grid** (one pseudobulk column per replicate×assay). Everything lives under per-assay `bb_dir/{assay_type}/`; the shared grid and sample sheet are duplicated into each sub-dir:
   - `{assay_type}/bb.tsv.gz` — the one shared bin grid (joint across assays; identical copy in each sub-dir).
   - `{assay_type}/sample_ids.tsv` — sample metadata (one row per replicate×assay; identical copy in each sub-dir).
   - `{assay_type}/bb.{Tallele,Aallele,Ballele}.npz` — per-assay allele count matrices (bins × cells).
   - `{assay_type}/bb.Xcount.npz` — per-assay native-count matrix per bb bin (same shape/order as the allele matrices). scATAC: deduped fragment counts from raw fragments. scRNA/VISIUM: UMI counts from the `process_rna_anndata` h5ad (gene → largest-overlap bin).
   - `{assay_type}/multi_snp.tsv.gz`, `{assay_type}/multi_snp.{Tallele,Aallele,Ballele}.npz` — per-assay multi-SNP groups.
   - `{assay_type}/barcodes.tsv.gz`, `{assay_type}/barcodes.full.tsv.gz` — per-assay cell barcodes.

To reuse genotyped and phased SNPs from bulk data, set `het_snp_vcf` to a prior run's `phase/phased_het_snps.vcf.gz`.

---

## Mode 3: `copytyping_preprocess`

Skips genotyping. Aggregates per-cell allele counts onto pre-computed BB blocks.

1. Sample sheet with non-bulk samples.
2. Set `workflow_mode: copytyping_preprocess`. Provide `het_snp_vcf` and `bb_file` in config.
3. Outputs in `bb_dir/{assay_type}/`:
   - `cnv_segments.tsv` — BB block annotations.
   - `bb.{Xcount,Tallele,Aallele,Ballele}.npz` — per-block count matrices. scATAC `Xcount` comes from raw 10x fragments; RNA-family from the `process_rna_anndata` h5ad.
   - `barcodes.tsv.gz`, `barcodes.full.tsv.gz` — cell barcode list.
   - `sample_ids.tsv` — sample metadata.

---

## Phasing

Set `phaser` in config to one of:

| Phaser | Config keys needed |
|--------|--------------------|
| `eagle` | `gmap_path` — full path to Eagle2 single gmap file (e.g., `genetic_map_hg38_withX.txt.gz`) |
| `shapeit` | `gmap_path` — full path with `{chrname}` placeholder (e.g., `.../chr{chrname}.b38.gmap.gz`) |
| `longphase` | `params_longphase` — `min_mapq`, `extra_params` (e.g., `"--pb"` for PacBio) |

Eagle and shapeit also require a phasing reference panel (`phasing_panel` in config). Longphase doesn't.

See [resources/README.md](../resources/README.md) for per-reference gmap sources (hg38 / chm13v2 / mm10).

Supported `reference_version` values: `hg19`, `hg38`, `chm13v2`, `mm10`.
