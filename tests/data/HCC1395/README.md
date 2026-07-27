# Integration test: HCC1395 / HCC1395BL (WGS + WES + ONT + HiFi, hg38, chr16, streaming)

Remote-only (no GitHub CI): even chr16-only, the eight BAMs stream from public servers.

## 1. Dataset

Real triple-negative breast cancer tumor-normal pair for `bulk_genotyping` on **hg38**,
**chr16 only**, in one multi-sample run across four assay types (all open-access, GRCh38,
`chr`-prefixed):

| dataset_id | assay_type | sample | platform | source |
|---|---|---|---|---|
| `N_ILL` / `T_ILL` | `bulkWGS` | normal / tumor | Illumina SR | SEQC2 Somatic Mutation WG, BWA-MEM to GRCh38.d1.vd1 |
| `N_WES` / `T_WES` | `bulkWES` | normal / tumor | Illumina SR | SEQC2 Somatic Mutation WG (WES/), same alignment |
| `N_ONT` / `T_ONT` | `bulkWGS-lr` | normal / tumor | ONT | CASTLE / Park 2024, minimap2 GRCh38 |
| `N_HIFI` / `T_HIFI` | `bulkWGS-lr` | normal / tumor | PacBio HiFi | CASTLE / Park 2024, pbmm2 GRCh38 |

All eight BAMs are read **directly by URL** (`remote_mode: stream`); the co-located index
(`##idx##`) lets htslib fetch only chr16, never the whole-genome/whole-exome BAMs. Nothing is
downloaded into the repo, and references stage under the gitignored `.test-run/`, so the
working tree stays clean.

One shared grid over all bulk assays: het SNPs are genotyped from the Illumina WGS normal
(`genotype_dataset_ids: ["N_ILL"]`, bcftools), phased with eagle against the 1000G n=3202
panel, then every dataset (WGS/WES/ONT/HiFi, tumor+normal) is piled up and depth-counted on
that grid. Each tumor's RDR divides by its matched-platform normal (`rdr_base_dataset_id`).

- SR WGS + WES: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/{WGS,WES}/
  (index is Picard-style `basename.bai`; SEQC2 BioProject PRJNA489865 / SRP162370)
- ONT + HiFi: https://storage.googleapis.com/brain-genomics-public/publications/park2024_deepsomatic/bams/
  (index is `basename.bam.bai`; CASTLE BioProject PRJNA1086849)

> [!NOTE]
> chr16 carries a clear within-chromosome CN/BAF transition: the q-arm is entirely lost with
> LOH in HCC1395 (16p retained). This is a strong *clonal* signal (good for CI); it is not a
> subclonal region. Masood et al., Genome Biology 25:163 (2024),
> doi:10.1186/s13059-024-03294-8 (6p / 16q / X entirely lost); Fang et al. 2021,
> doi:10.1038/s41587-021-00993-6. `bulkWES` depth carries capture-enrichment structure, so its
> post-correction RDR is noisier; inspect `qc/rd_correction.bulkWES.pdf`.

## 2. Remote server (conda)

Needs an env with `snakemake>=9` + `mamba`, and `bcftools`/`samtools`/`tabix`/`bgzip`/`wget`
on PATH for `prepare_refs.sh`. `run_pipeline.sh` builds per-rule envs itself (`--use-conda`).
Set `local-storage-prefix` (profile) and `RUN` to a scratch disk. In stream mode only chr16
byte ranges transfer, so peak disk stays small.

```sh
bash tests/data/HCC1395/prepare_refs.sh
RUN=/scratch/HCC1395 CORES=<N> bash tests/data/HCC1395/run_pipeline.sh
RUN=/scratch/HCC1395 bash tests/data/HCC1395/check_outputs.sh
```

To run genome-wide (autosomes 1-22) instead of chr16, stage whole-genome references
(`CHROMS`) and use the whole-genome runner; `check_outputs.sh` takes the matching `RUN`:

```sh
CHROMS="$(seq 1 22)" bash tests/data/HCC1395/prepare_refs.sh
RUN=/scratch/HCC1395_wg CORES=<N> bash tests/data/HCC1395/run_pipeline.whole_genome.sh
RUN=/scratch/HCC1395_wg bash tests/data/HCC1395/check_outputs.sh
```
