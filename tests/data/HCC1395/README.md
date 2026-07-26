# Integration test: HCC1395 / HCC1395BL (SEQC2, hg38 Illumina)

## 1. Dataset

Real triple-negative breast cancer tumor-normal WGS pair for `bulk_genotyping` on **hg38**.
SEQC2 Illumina NovaSeq, BWA-MEM to GRCh38.d1.vd1 (chr-prefixed, no ALT contigs); whole-genome
BAMs (~114/127 GiB) stream from NCBI ReferenceSamples by URL and are region-subset via the
remote index (`subset_remote_alignment`), so the full BAM never downloads. Genotyped with
bcftools, phased with eagle against the 1000G n=3202 panel. Uses **chr6 + chr16** (6p and 16q
entire-arm losses with LOH) -- clear within-chromosome CN/BAF transitions.

- Directory + md5s: https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WGS/
- SEQC2: BioProject PRJNA489865 / SRA SRP162370.
- CNAs: Gong, Fang et al., Genome Biology 25:158 (2024), doi:10.1186/s13059-024-03294-8

Index is Picard-style `basename.bai` (not `basename.bam.bai`); `samples.json` points
`alignment_index` at the exact `.bai` URL.

## 2. GitHub CI (or `act` + Docker)

```sh
gh workflow run it-HCC1395.yml   # GitHub Actions (workflow_dispatch)

# locally in Docker via act (.actrc pins catthehacker/ubuntu:act-24.04 + linux/amd64; ~16 GB VM for eagle)
act workflow_dispatch -W .github/workflows/it-HCC1395.yml
```

## 3. Remote server (conda, no Docker/act)

Needs an env with `snakemake>=9` + `mamba`, and `bcftools`/`samtools`/`tabix`/`bgzip`/`wget`
on PATH for `prepare_refs.sh`. `run_pipeline.sh` builds per-rule envs itself (`--use-conda`).
Runtime (staged refs + outputs) lands under the gitignored `.test-run/HCC1395/`.

```sh
bash tests/data/HCC1395/prepare_refs.sh
CORES=<N> bash tests/data/HCC1395/run_pipeline.sh
bash tests/data/HCC1395/check_outputs.sh
```
