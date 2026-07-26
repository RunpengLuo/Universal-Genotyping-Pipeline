# Integration test: HCC1395 / HCC1395BL (SEQC2, hg38 Illumina)

Remote-only (no GitHub CI): the whole-genome BAMs are too large for a hosted runner.

## 1. Dataset

Real triple-negative breast cancer tumor-normal WGS pair for `bulk_genotyping` on **hg38**.
SEQC2 Illumina NovaSeq, BWA-MEM to GRCh38.d1.vd1 (chr-prefixed, no ALT contigs). Whole-genome
BAMs (~114/127 GiB) stream from NCBI ReferenceSamples by URL and are fetched in full by
Snakemake `storage()` into `local-storage-prefix` (point it at scratch). Genotyped with
bcftools, phased with eagle against the 1000G n=3202 panel. Uses **chr6 + chr16** (6p and 16q
entire-arm losses with LOH) -- clear within-chromosome CN/BAF transitions.

- Directory + md5s: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/
- SEQC2: BioProject PRJNA489865 / SRA SRP162370.
- CNAs: Gong, Fang et al., Genome Biology 25:158 (2024), doi:10.1186/s13059-024-03294-8

Index is Picard-style `basename.bai` (not `basename.bam.bai`); `samples.json` points
`alignment_index` at the exact `.bai` URL.

## 2. Remote server (conda)

Needs an env with `snakemake>=9` + `mamba`, and `bcftools`/`samtools`/`tabix`/`bgzip`/`wget`
on PATH for `prepare_refs.sh`. `run_pipeline.sh` builds per-rule envs itself (`--use-conda`).
Set `local-storage-prefix` (profile) and `RUN` to a large scratch disk -- storage() downloads
the whole BAMs. Runtime lands under `$RUN`.

```sh
bash tests/data/HCC1395/prepare_refs.sh
RUN=/scratch/HCC1395 CORES=<N> bash tests/data/HCC1395/run_pipeline.sh
RUN=/scratch/HCC1395 bash tests/data/HCC1395/check_outputs.sh
```
