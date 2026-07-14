# Tutorial: Bulk WGS Genotyping (hg38)

`bulk_genotyping` on paired normal/tumor bulk WGS aligned to GRCh38. Install first: [README](../../README.md#installation).

## 1. Sample file

Schema: [sample_sheet.md](../sample_sheet.md).

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

Validate: [README](../../README.md#quick-start).

## 2. Config

Copy `resources/templates/config.yaml`. Every key: [reference.md](../reference.md#config-keys).

```yaml
workflow_mode: "bulk_genotyping"
assay_types: ["bulkWGS"]

sample_id: HT001
sample_file: /path/to/samples.json
chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

reference_version: hg38
reference: /path/to/reference.fasta
genome_size: resources/data/hg38.regions.bed
region_bed: resources/data/hg38.regions.bed
window_bed: resources/data/windows.1kbp.hg38.bed.gz
blacklist_bed: resources/data/hg38-blacklist.v2.bed.gz
gtf_file: /path/to/gencode.v38.annotation.gtf.gz

params_combine_counts:
  min_snp_reads: [500, 1000, 2000, 3000]
```

### Reference resources

| Key | Resource |
|-----|----------|
| `window_bed` | `windows.1kbp.hg38.bed.gz` — 1 kb windows with GC / mappability / replication timing (in `resources/data/`) |
| `blacklist_bed` | `hg38-blacklist.v2.bed.gz` — ENCODE blacklist v2 (in `resources/data/`) |
| `region_bed`, `genome_size` | `hg38.regions.bed` — chromosome regions, minus centromeres and acrocentric p-arms (13/14/15/21/22) (in `resources/data/`) |
| `gtf_file` | GENCODE v38 ([download](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz)) |

### SNP and phasing panels

[1kGP n=3,202 high-coverage](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/), prepared with bcftools / bgzip / tabix:

```bash
bash resources/scripts/process_1kGP_3202_panel.sh --ref hg38 /path/to/1kGP_3202
```

```yaml
snp_panel: /path/to/snps.vcf.gz
snp_targets: /path/to/target_positions
phasing_panel: /path/to/phasing_panel
```

### Phasing

SHAPEIT5 ([GitHub](https://github.com/odelaneau/shapeit5)) is recommended; Eagle2 ([download](https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz)) also works. Each ships its own genetic map.

```yaml
phaser: "shapeit"
gmap_path: "/path/to/shapeit5/resources/maps/b38/chr{chrname}.b38.gmap.gz"
```

```yaml
phaser: "eagle"
gmap_path: "/path/to/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
```

## 3. Run

Set `cores` in `profile/config.yaml`, then:

```bash
snakemake --profile /path/to/profile/ \
    -s /path/to/workflow/Snakefile \
    --configfile /path/to/my_config.yaml \
    --directory /path/to/output_dir \
    --config sample_file=/path/to/samples.json sample_id=HT001
```

## 4. Results

`<out_dir>/<bb_dir>/MSR{msr}/{stream}/` (here `stream` = `bulkWGS`): `bb.tsv.gz`, `bb.{Tallele,Aallele,Ballele,depth,rdr}.npz`, `sample_ids.tsv`. Columns: [reference.md](../reference.md#outputs).

QC: `<qc_dir>/rd_correction.bulkWGS.pdf` (bias correction), `combine_counts.bulkWGS.MSR{msr}.pdf` (one per `min_snp_reads`; compare BAF/RDR plots to pick a bin size).
