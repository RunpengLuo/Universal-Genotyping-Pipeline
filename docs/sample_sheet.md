# Sample File Format

The sample file (`sample_file` in config) is JSON: a list of records, one per
**(`dataset_id`, `assay_type`)**. Only records whose `sample_id` matches the config `sample_id` and
whose `assay_type` is in `assay_types` are used, so one file can describe many individuals and assays.

Optional keys are omitted, never blanked. Unknown keys are rejected, so typos fail loudly.

`key?` marks an optional key; `|` separates alternatives; `<path|url>` is a local path or an
`http(s)` URL.

```
{
  "version": 1,
  "samples": [
    {
      "sample_id":             <str>,
      "dataset_id":            <str>,
      "rdr_base_dataset_id"?:  <str>,
      "passage_id"?:           <str>,
      "assay_type":            "bulkWGS" | "bulkWGS-lr" | "bulkWES" |
                               "scRNA" | "scATAC" | "VISIUM" | "VISIUM3prime",
      "sample_type":           "normal" | "tumor",
      "platform"?:             <str>,
      "reference_version"?:    <str>,
      "files": {
        "alignment":         <path|url>,
        "alignment_index":   <path|url>,
        "barcodes"?:         <path|url>,
        "fragments"?:        <path|url>,
        "matrix_h5"?:        <path|url>,
        "tissue_positions"?: <path|url>,
        "scalefactors"?:     <path|url>,
        "image_hires"?:      <path|url>,
        "image_lowres"?:     <path|url>
      },
      "meta"?:                 { <str>: <any>, ... }
    },
    ...
  ]
}
```

A key optional here may still be required by the record's `assay_type` — see the tables below. A
10x multiome dataset is two records sharing one `dataset_id`, one `scRNA` and one `scATAC`.

## Record keys

| Key | Type | Required | Description |
|-----|------|----------|-------------|
| `sample_id` | string | Yes | Individual ID; matched against `sample_id` in config. |
| `dataset_id` | string | Yes | Dataset (replicate) ID; unique per `(dataset_id, assay_type)`. A multiome pair shares one. |
| `rdr_base_dataset_id` | string | No | Bulk tumor only: the `dataset_id` this tumor's depth is divided by for RDR. Omit -> median normalization. |
| `passage_id` | string | No | Cell-line passage (e.g. `p23`). |
| `assay_type` | string | Yes | `bulkWGS`, `bulkWGS-lr`, `bulkWES`, `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime`. |
| `sample_type` | string | Yes | `normal` or `tumor`. |
| `platform` | string | No | Sequencing platform (e.g. `Illumina PCR-free`, `PacBio Revio HiFi`). Free text; a tumor and its `rdr_base_dataset_id` should agree, so RDR cancels platform GC/mappability bias. |
| `reference_version` | string | No | Reference the `bam` is aligned to (e.g. `GRCh38-GIABv3`). Free text and independent of the config key of the same name; records selected for one run must be coordinate-compatible, see [Multi-dataset samples](#multi-dataset-samples). |
| `cancer_type` | string | No | Tumor type of the individual (e.g. `BRCA_Basal`). Free text. |
| `files` | object | Yes | Input files; see below. |

Only `sample_id`, `dataset_id`, `rdr_base_dataset_id`, `assay_type`, `sample_type`, and `files`
change what the pipeline does. `passage_id`, `platform`, `reference_version`, and `cancer_type` are
provenance: coerced to string, but no rule reads them.

A record may carry **any other key** — the table above is not a whitelist. Unrecognized keys are kept
verbatim and never read, so a record can hold provenance freely (e.g. `source_url`, `clone_id`,
`coverage`, or a nested `meta` object of your own).

## Files

Every file is named explicitly — no directory is ever an input, so any of them may be remote. Each
value becomes a tracked Snakemake input.

| Key | Required for | Typical 10x source |
|-----|--------------|--------------------|
| `alignment` | every assay type | `outs/*_possorted_bam.bam` |
| `alignment_index` | every assay type | `outs/*_possorted_bam.bam.bai` |
| `barcodes` | `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime` | `outs/filtered_feature_bc_matrix/barcodes.tsv.gz` |
| `fragments` | `scATAC` | `outs/atac_fragments.tsv.gz` |
| `matrix_h5` | `scRNA`, `VISIUM`, `VISIUM3prime` | `outs/filtered_feature_bc_matrix.h5` |
| `tissue_positions` | `VISIUM`, `VISIUM3prime` | `outs/spatial/tissue_positions.csv` (or `tissue_positions_list.csv`) |
| `scalefactors` | `VISIUM`, `VISIUM3prime` | `outs/spatial/scalefactors_json.json` |
| `image_hires` | `VISIUM` | `outs/spatial/tissue_hires_image.png` |
| `image_lowres` | `VISIUM` | `outs/spatial/tissue_lowres_image.png` |

`alignment` is a `.bam` or a `.cram`; `alignment_index` is its `.bai` or `.crai`. Both are always
required — nothing is inferred from a filename. The remaining requirements apply only to
`single_cell_genotyping` and `copytyping_preprocess` runs. Supporting new, non-10x data means adding
a key here, not a column to every record.

The table is the set of files the pipeline **reads**, not a whitelist of what may appear. Any other
key in `files` — an unlisted name such as `fastq_r1`, or a listed one the assay does not consume — is
dropped at parse time with a `NOTE:` on stderr: no rule inputs it, it never claims a download slot,
and `--check-files` does not stat it. Record-only paths are therefore free to live in `files`. The
flip side: a typo'd key (`alignment_idx`) is silently dropped and then reported as the required key
being missing.

Spatial records are read by squidpy: the pipeline stages the named files into a Space Ranger layout
and calls `squidpy.read.visium()` on it, so `obsm["spatial"]` and `uns["spatial"]` (scalefactors and
both tissue images) are populated exactly as they would be from a local `outs/`. `VISIUM3prime`
carries no images — squidpy cannot load them for 3' data — so `image_hires` / `image_lowres` are
ignored on those records.

## Remote inputs

Every `files` value accepts an `http(s)` URL in place of a local path, detected per value — a record
may mix remote and local files. Snakemake fetches each one once into `.snakemake/storage/`, shares it
across every rule that needs it, and deletes it once no job needs it. Requires
`snakemake-storage-plugin-http`, in `environment.yaml`. There is no streaming mode: a fetched file is
held whole on disk.

Retrieval is throttled by the `downloads` resource, set in `profile/config.yaml`:

| Setting | Default | Effect |
|---------|---------|--------|
| `resources: downloads=N` | 2 | At most N jobs that must fetch a URL run at once. Peak local disk is about N x the largest alignment, and the host sees at most N concurrent requests. |
| `keep-storage-local-copies` | false | Keep fetched files instead of deleting them after the last job that needs them. |
| `local-storage-prefix` | `.snakemake/storage` | Where fetched files land; point at a scratch disk if the working dir is small. |

Raise `downloads` for a fast link and a large disk; set it to 1 to serialize retrieval. Public hosts
may rate-limit: the NCBI GIAB host answers **503** to a dozen concurrent requests, so a high value
will fail rather than go faster.

`alignment` and `alignment_index` are fetched together, so an index published under any name works.

## Multi-dataset samples

One `sample_id` may hold many datasets of the same individual. Genotyping and phasing run once per
`sample_id` and are shared by every record, which is correct only when:

- **Records share coordinates.** `reference_version` values may be mixed only if their primary chromosomes
  are identical: `GRCh38` and `GRCh38-GIABv3` qualify (GIABv3 masks GRC false duplications and
  appends decoys, leaving chr1-22/X/Y unchanged). Never mix GRCh38 with GRCh37 or CHM13 — use one
  run per incompatible build.
- **`dataset_id` disambiguates origin.** A tumor and its matched normal often share a source dataset
  id, and dataset ids repeat across builds; encode both, e.g. `Tp23-Element-2-GIABv3`.
- **Each tumor's `rdr_base_dataset_id` names a normal of the same `platform`** (and, where possible,
  the same `reference_version`), so RDR cancels platform GC and mappability bias.

For GIAB HG008, `~/Research/datasets/GIAB/giab_samplesheet.py` generates this file from the dataset manifests.

## Legacy TSV

A `.tsv` sample sheet is still accepted, with the original columns: `SAMPLE`, `REP_ID`,
`assay_type`, `sample_type`, `PATH_to_bam` (required), plus `RDR_BASE_REP_ID`, `passage`,
`PATH_to_barcodes`, `PATH_to_10x_ranger`. `parse_sample_file()` dispatches on the extension and maps
the columns onto the record schema above.

TSV has no column for the index, so `alignment_index` is derived as `PATH_to_bam` + `.bai` / `.crai`.
The single-cell files likewise have no columns, so they are re-derived from `PATH_to_10x_ranger` using Cell/Space Ranger's default layout: `filtered_feature_bc_matrix.h5`,
`atac_fragments.tsv.gz` (else `fragments.tsv.gz`), and `spatial/tissue_positions.csv` (else
`tissue_positions_list.csv`), `spatial/scalefactors_json.json`, `spatial/tissue_hires_image.png`,
`spatial/tissue_lowres_image.png`. Those names live in `config/const.py` (`RANGER_*`), so a Ranger
release that renames a file is a one-file change. Consequences:

- `PATH_to_10x_ranger` must be a **local** directory: a directory cannot be fetched, so single-cell
  TSV sheets cannot use remote inputs. `PATH_to_bam` may still be a URL.
- A non-default filename inside `outs/` cannot be expressed. Use JSON.

New sheets should use JSON; TSV exists so existing sheets keep running.

## Validation

`resources/scripts/validate_sample_file.py <sample_file>` runs these checks without starting the
workflow (`--check-files` also stats every local path; `--sample-id` / `--workflow-mode` narrow the
scope). The same checks run before the DAG is built; errors name the record by `sample_id`,
`dataset_id`, and `assay_type`:

- Missing required key; `assay_type` / `sample_type` out of range.
- A required file missing for the assay type. (Extra record keys and extra `files` keys are not
  errors: they are kept and ignored, respectively.)
- Duplicate `(dataset_id, assay_type)`; a `dataset_id` with >2 assays or a 2-assay `dataset_id` that
  is not an `scRNA` + `scATAC` pair; duplicate `dataset_id` in bulk mode.
- `rdr_base_dataset_id` on a non-tumor record, naming itself, or naming a `dataset_id` absent from
  the `sample_id`.
