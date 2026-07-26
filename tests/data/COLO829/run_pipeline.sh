#!/usr/bin/env bash
# Run a tests/data/<case> end-to-end to bb_dir. A step of the CI job (run via act locally).
#
# Runpeng Luo (2026-07-24)
#
# Dependencies: snakemake (>=9) on PATH, --use-conda per-rule envs (via profile/).
#   Run inside the genotyping-env env (the CI/act job activates it).
#
# Usage:
#   [RUN=<out_dir>] [CONFIG=<config>] [CORES=<n>] bash run_pipeline.sh
#           RUN     # output root (default: <repo>/.test-run/<case>)
#           CONFIG  # config file (default: <case>/config.yaml)
#           CORES   # snakemake --cores (default: 4)
# Notes:
#   cwd is forced to the repo root so the config's repo-relative reference paths
#   resolve; output dirs are redirected under RUN. CASE is inferred from this
#   script's location, so the file is identical across cases.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
cd "$REPO"

CASE="${HERE#"$REPO"/}"
CASENAME="$(basename "$HERE")"
RUN="${RUN:-$REPO/.test-run/$CASENAME}"
CONFIG="${CONFIG:-$CASE/config.yaml}"
CORES="${CORES:-4}"

snakemake -s workflow/Snakefile \
  --configfile "$CONFIG" \
  --profile profile/ --cores "$CORES" \
  --config \
    snp_dir="$RUN/snps" phase_dir="$RUN/phase" pileup_dir="$RUN/pileup" \
    allele_dir="$RUN/allele" bb_dir="$RUN/bb" qc_dir="$RUN/qc" \
    log_dir="$RUN/logs" aux_dir="$RUN/aux" bench_dir="$RUN/benchmarks"
