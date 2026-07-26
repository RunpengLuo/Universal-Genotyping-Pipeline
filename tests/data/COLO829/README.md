# Integration test: COLO829 / COLO829BL (hg38, ONT long-read)

Remote-only (no GitHub CI): the whole-genome BAMs are too large for a hosted runner.

## 1. Dataset

Real melanoma tumor-normal WGS pair for `bulk_genotyping` on **hg38**, long-read (Oxford
Nanopore). ONT Open Data `colo829_2024.03`: Dorado sup basecalls, minimap2 to GRCh38
(chr-prefixed). Whole-genome BAMs (~135/167 GiB) stream from the AWS `ont-open-data` bucket
by URL and are fetched in full by Snakemake `storage()` into `local-storage-prefix` (point it
at scratch). Long-read path: `bulkWGS-lr` + longphase (`--ont`), het SNPs at 1000G n=3202
targets. Default config uses **chr8 + chr10** (trisomy 8; chr10 loss/LOH + PTEN del) -- clear
allelic imbalance. Variant `config.full.yaml` (+ `prepare_refs.full.sh`) runs all 22 autosomes.

- Release notes: https://labs.epi2me.io/colo-2024.03/
- CNAs: Velazquez-Villarreal et al., Commun Biol 3:318 (2020), doi:10.1038/s42003-020-1044-8

## 2. Remote server (conda)

Needs an env with `snakemake>=9` + `mamba`, and `bcftools`/`samtools`/`tabix`/`bgzip`/`wget`
on PATH for `prepare_refs.sh`. `run_pipeline.sh` builds per-rule envs itself (`--use-conda`).
Set `local-storage-prefix` (profile) and `RUN` to a large scratch disk -- storage() downloads
the whole BAMs. Runtime lands under `$RUN`.

```sh
# chr8 + chr10
bash tests/data/COLO829/prepare_refs.sh
RUN=/scratch/COLO829 CORES=<N> bash tests/data/COLO829/run_pipeline.sh
RUN=/scratch/COLO829 bash tests/data/COLO829/check_outputs.sh

# full autosomes (chr1-chr22)
bash tests/data/COLO829/prepare_refs.full.sh
RUN=/scratch/COLO829_full CONFIG=tests/data/COLO829/config.full.yaml CORES=<N> \
  bash tests/data/COLO829/run_pipeline.sh
RUN=/scratch/COLO829_full bash tests/data/COLO829/check_outputs.sh
```
