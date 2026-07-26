# Integration test: hatchet2_chr22_simulation

## 1. Dataset

HATCHet demo-complete simulation: 1 matched normal + 3 tumors, **chr22, hg19/GRCh37**,
`bulk_genotyping`. Alignments stream from Zenodo by URL (no local copy). Default config
genotypes with bcftools and phases with eagle against the 1000GP Phase3 chr22 panel;
uses the pre-built hg19 window BED (GC+MAP+REPLI baked in) + blacklist + GENCODE v19.

- Alignments: https://zenodo.org/records/4046906 -- Panel: https://zenodo.org/records/6709541
- https://raphael-group.github.io/hatchet/examples/demo-complete/demo-complete.html

Variant `config.dbsnp151.yaml` (+ `prepare_refs.dbsnp151.sh`): same reads genotyped at
dbSNP151 chr22 positions and phased with shapeit5.

## 2. GitHub CI (or `act` + Docker)

```sh
gh workflow run it-hatchet2_chr22_simulation.yml           # GitHub Actions (workflow_dispatch)
gh workflow run it-hatchet2_chr22_simulation_dbsnp151.yml  # dbSNP151 + shapeit5 variant

# locally in Docker via act (.actrc pins catthehacker/ubuntu:act-24.04 + linux/amd64; ~16 GB VM)
act workflow_dispatch -W .github/workflows/it-hatchet2_chr22_simulation.yml
```

## 3. Remote server (conda, no Docker/act)

Needs an env with `snakemake>=9` + `mamba`, and `bcftools`/`samtools`/`tabix`/`bgzip`/`wget`
on PATH for `prepare_refs.sh`. `run_pipeline.sh` builds per-rule envs itself (`--use-conda`).
Runtime (staged refs + outputs) lands under the gitignored `.test-run/hatchet2_chr22_simulation/`.

```sh
bash tests/data/hatchet2_chr22_simulation/prepare_refs.sh
CORES=<N> bash tests/data/hatchet2_chr22_simulation/run_pipeline.sh
bash tests/data/hatchet2_chr22_simulation/check_outputs.sh

# dbSNP151 + shapeit5 variant
bash tests/data/hatchet2_chr22_simulation/prepare_refs.dbsnp151.sh
RUN="$PWD/.test-run-dbsnp151/hatchet2_chr22_simulation" \
  CONFIG=tests/data/hatchet2_chr22_simulation/config.dbsnp151.yaml CORES=<N> \
  bash tests/data/hatchet2_chr22_simulation/run_pipeline.sh
RUN="$PWD/.test-run-dbsnp151/hatchet2_chr22_simulation" \
  bash tests/data/hatchet2_chr22_simulation/check_outputs.sh
```
