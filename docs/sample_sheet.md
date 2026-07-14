# Sample File Format

Sample file is a JSON-format configuration file contains dataset records over multiple patients and patient-specific datasets. Template can be found at here: [templates](../resources/templates/).


## Schema

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
      "cancer_type"?:          <str>,
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
      }
    },
    ...
  ]
}
```

## Fields

| Key | Required | Description |
|-----|----------|-------------|
| `sample_id` | Yes | Patient ID |
| `dataset_id` | Yes | Dataset ID, unique with each patient |
| `assay_type` | Yes | `bulkWGS` \| `bulkWGS-lr` \| `bulkWES` \| `scRNA` \| `scATAC` \| `VISIUM` \| `VISIUM3prime`. |
| `sample_type` | Yes | `normal` \| `tumor`. |
| `files` | Yes | Input files; see [Files](#files). |
| `rdr_base_dataset_id` | No | Normal `dataset_id` used as the RDR baseline; see [RDR normalization](workflow.md#bulk_genotyping). |
| `passage_id` | No | e.g. `p23`. |
| `platform` | No | e.g. `Illumina PCR-free`. |
| `reference_version` | No | e.g. `GRCh38-GIABv3`; see [Workflow](workflow.md). |
| `cancer_type` | No | e.g. `PDAC`. |

## Files

| Key | Required for | Example |
|-----|--------------|--------------------|
| `alignment` | every assay type | `outs/*_possorted_bam.bam` |
| `alignment_index` | every assay type | `outs/*_possorted_bam.bam.bai` |
| `barcodes` | `scRNA`, `scATAC`, `VISIUM`, `VISIUM3prime` | `outs/filtered_feature_bc_matrix/barcodes.tsv.gz` |
| `fragments` | `scATAC` | `outs/atac_fragments.tsv.gz` |
| `matrix_h5` | `scRNA`, `VISIUM`, `VISIUM3prime` | `outs/filtered_feature_bc_matrix.h5` |
| `tissue_positions` | `VISIUM`, `VISIUM3prime` | `outs/spatial/tissue_positions.csv` |
| `scalefactors` | `VISIUM`, `VISIUM3prime` | `outs/spatial/scalefactors_json.json` |
| `image_hires` | `VISIUM` | `outs/spatial/tissue_hires_image.png` |
| `image_lowres` | `VISIUM` | `outs/spatial/tissue_lowres_image.png` |

## Features
### Support multiome dataset
A 10x Epi Multiome dataset is two records sharing one `dataset_id`, one `scRNA` and one `scATAC`.

### Support remote files
A input file could be a local disk path or an `http(s)` URL to a FTP server. Note that `snakemake-storage-plugin-http` (in `environment.yaml`) is required. Remote files will be downloaded to `local-storage-prefix` when request by a rule and deleted afterwards. See `profile/config.yaml` for more information.

### Validate sample file
Use `resources/scripts/validate_sample_file.py <sample_file>` to valid the sample file's format.

## Legacy TSV
Still accepted; `parse_sample_file()` dispatches on the extension and maps the columns onto the schema above. Columns: `SAMPLE`, `REP_ID`, `assay_type`, `sample_type`, `PATH_to_bam` (required), plus `RDR_BASE_REP_ID`, `passage`, `PATH_to_barcodes`, `PATH_to_10x_ranger`.

TSV has no index column, so `alignment_index` is derived as `PATH_to_bam` + `.bai` / `.crai`. Single-cell files have no columns either and are re-derived from `PATH_to_10x_ranger` using Cell/Space Ranger's default layout (names live in `config/const.py` as `RANGER_*`, so a Ranger rename is a one-file change). Consequences:

- `PATH_to_10x_ranger` must be a **local** directory — a directory cannot be fetched, so single-cell TSV sheets cannot use remote inputs. `PATH_to_bam` may still be a URL.
- A non-default filename inside `outs/` cannot be expressed.
