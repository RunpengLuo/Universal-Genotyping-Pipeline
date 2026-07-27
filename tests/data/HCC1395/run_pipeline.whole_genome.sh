#!/usr/bin/env bash
# Run the HCC1395 case genome-wide (autosomes 1-22) instead of the chr16-only default.
#
# Runpeng Luo (2026-07-27)
#
# Dependencies: snakemake (>=9) on PATH, --use-conda per-rule envs (via profile/).
#   Run inside the genotyping-env env. Stage whole-genome references first:
#     CHROMS="$(seq 1 22)" bash tests/data/HCC1395/prepare_refs.sh
#
# Usage:
#   [RUN=<out_dir>] [CONFIG=<config>] [CORES=<n>] [PROFILE=<dir>] bash run_pipeline.whole_genome.sh
#           RUN      # output root (default: <repo>/.test-run/HCC1395_wg)
#           CONFIG   # config file (default: <case>/config.yaml)
#           CORES    # snakemake --cores (default: 4)
#           PROFILE  # snakemake --profile dir (default: profile/, repo-relative)
# Notes:
#   Same config and sample file as run_pipeline.sh; only `chromosomes` is overridden
#   to the 22 autosomes. All eight BAMs stream by URL (remote_mode: stream), so no
#   whole-file download - only the per-chromosome byte ranges transfer.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

CASE="${HERE#"$REPO"/}"
RUN="${RUN:-$REPO/.test-run/HCC1395_wg}"
CONFIG="${CONFIG:-$CASE/config.yaml}"
CORES="${CORES:-4}"
PROFILE="${PROFILE:-profile/}"

CHROMS="$(seq -s , 1 22)"

snakemake -s workflow/Snakefile \
  --configfile "$CONFIG" \
  --profile "$PROFILE" --cores "$CORES" \
  --config \
    chromosomes="[${CHROMS}]" \
    snp_dir="$RUN/snps" phase_dir="$RUN/phase" pileup_dir="$RUN/pileup" \
    allele_dir="$RUN/allele" bb_dir="$RUN/bb" qc_dir="$RUN/qc" \
    log_dir="$RUN/logs" aux_dir="$RUN/aux" bench_dir="$RUN/benchmarks"
