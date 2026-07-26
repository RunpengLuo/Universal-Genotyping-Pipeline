# Integration test: COLO829 / COLO829BL (hg38, ONT long-read)

## 1. Dataset

Real melanoma tumor-normal WGS pair for `bulk_genotyping` on **hg38**, long-read (Oxford
Nanopore). ONT Open Data `colo829_2024.03`: Dorado sup basecalls, minimap2 to GRCh38
(chr-prefixed); whole-genome BAMs stream from the AWS `ont-open-data` bucket by URL and are
region-subset via the remote index (`subset_remote_alignment`), so the full BAM never
downloads. Long-read path: `bulkWGS-lr` + longphase (`--ont`), het SNPs at 1000G n=3202
targets. Default config uses **chr8 + chr10** (trisomy 8; chr10 loss/LOH + PTEN del) --
clear allelic imbalance.

- Release notes: https://labs.epi2me.io/colo-2024.03/
- CNAs: Velazquez-Villarreal et al., Commun Biol 3:318 (2020), doi:10.1038/s42003-020-1044-8

Variant `config.full.yaml` (+ `prepare_refs.full.sh`): all 22 autosomes (full COLO829 CNA
landscape). **No CI** -- subsetting every autosome fetches ~the whole ONT BAM
(~130-167 GiB/sample), above a hosted runner's disk; remote-server only.

## 2. GitHub CI (or `act` + Docker)

Only the chr8+chr10 config runs in CI (the full variant is too large -- see section 3).

```sh
gh workflow run it-COLO829.yml   # GitHub Actions (workflow_dispatch)

# locally in Docker via act (.actrc pins catthehacker/ubuntu:act-24.04 + linux/amd64)
act workflow_dispatch -W .github/workflows/it-COLO829.yml
```

## 3. Remote server (conda, no Docker/act)

Needs an env with `snakemake>=9` + `mamba`, and `bcftools`/`samtools`/`tabix`/`bgzip`/`wget`
on PATH for `prepare_refs.sh`. `run_pipeline.sh` builds per-rule envs itself (`--use-conda`).
Runtime (staged refs + outputs) lands under the gitignored `.test-run/COLO829*/`.

```sh
bash tests/data/COLO829/prepare_refs.sh
CORES=<N> bash tests/data/COLO829/run_pipeline.sh
bash tests/data/COLO829/check_outputs.sh

# full-autosome variant (large disk; remote only)
bash tests/data/COLO829/prepare_refs.full.sh
RUN="$PWD/.test-run/COLO829_full" \
  CONFIG=tests/data/COLO829/config.full.yaml CORES=<N> \
  bash tests/data/COLO829/run_pipeline.sh
RUN="$PWD/.test-run/COLO829_full" bash tests/data/COLO829/check_outputs.sh
```
