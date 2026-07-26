#!/usr/bin/env bash
# Stage hg38 chr8+chr10 reference + het-SNP targets for the COLO829 (ONT) test case.
#
# Runpeng Luo (2026-07-26)
#
# Dependencies (on PATH): curl, wget, samtools, bcftools, bgzip, tabix (+ awk/gzip).
#   A minimal env suffices, e.g.
#   `mamba create -n genotyping-prep -c bioconda -c conda-forge bcftools samtools htslib wget`.
#
# Usage:
#   bash prepare_refs.sh [out_dir]
#           out_dir  # default: <repo>/.test-run/COLO829/reference (gitignored)
# Outputs (match tests/data/COLO829/config.yaml keys):
#   hg38.fa (+ .fai)                          -> reference (chr8 + chr10)
#   gencode.hg38.gtf.gz                        -> gtf_file (chr8 + chr10)
#   target_positions/target.chr{8,10}.pos.gz  -> snp_targets (1000G n=3202 het sites)
# Notes/References:
#   Reference: UCSC hg38 per-chromosome FASTA (GRCh38 primary; matches the chr-prefixed
#     GRCh38 ONT BAMs). Phasing is read-backed (longphase), so no phasing panel / genetic
#     map is staged; only het-SNP `snp_targets` come from the 1000G panel.
#   1000G n=3202 phased panel (GRCh38): resources/scripts/process_1kGP_3202_panel.sh.
#   GENCODE v44: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
#   Prebuilt genome-wide hg38 window BED + blacklist (resources/data/) are used directly.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../../.." && pwd)"
OUT="${1:-${REPO}/.test-run/COLO829/reference}"
mkdir -p "${OUT}"

UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
GENCODE_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

echo "[1/3] reference FASTA (hg38 chr8 + chr10)"
if [ ! -f "${OUT}/hg38.fa.fai" ]; then
  : > "${OUT}/hg38.fa"
  for c in chr8 chr10; do
    curl -fsSL "${UCSC}/${c}.fa.gz" | gunzip -c >> "${OUT}/hg38.fa"
  done
  samtools faidx "${OUT}/hg38.fa"
fi

echo "[2/3] GENCODE v44 GTF (chr8 + chr10)"
if [ ! -f "${OUT}/gencode.hg38.gtf.gz" ]; then
  curl -fsSL "${GENCODE_GTF}" | gunzip -c \
    | awk -F'\t' '$1=="chr8" || $1=="chr10"' \
    | gzip -c > "${OUT}/gencode.hg38.gtf.gz"
fi

echo "[3/3] 1000G n=3202 het-SNP targets (chr8, chr10)"
if [ ! -f "${OUT}/target_positions/target.chr10.pos.gz.tbi" ]; then
  bash "${REPO}/resources/scripts/process_1kGP_3202_panel.sh" \
    --ref hg38 --chroms "8 10" "${OUT}"
fi

echo "done: staged COLO829 reference under ${OUT}"
