#!/usr/bin/env bash
# Stage hg38 chr1-chr22 reference + het-SNP targets for the COLO829 full-autosome case.
#
# Runpeng Luo (2026-07-26)
#
# Dependencies (on PATH): curl, wget, samtools, bcftools, bgzip, tabix (+ awk/gzip).
#   A minimal env suffices, e.g.
#   `mamba create -n genotyping-prep -c bioconda -c conda-forge bcftools samtools htslib wget`.
#
# Usage:
#   bash prepare_refs.full.sh [out_dir]
#           out_dir  # default: <repo>/.test-run/COLO829_full/reference (gitignored)
# Outputs (match tests/data/COLO829/config.full.yaml keys):
#   hg38.fa (+ .fai)                              -> reference (chr1..chr22)
#   gencode.hg38.gtf.gz                           -> gtf_file (chr1..chr22)
#   target_positions/target.chr{1..22}.pos.gz    -> snp_targets (1000G n=3202 het sites)
# Notes/References:
#   Full-autosome variant of prepare_refs.sh (chr8+chr10). Same sources; genome-wide.
#   UCSC hg38 per-chromosome FASTA (GRCh38 primary; matches the chr-prefixed ONT BAMs).
#   Phasing is read-backed (longphase), so no phasing panel / genetic map is staged.
#   1000G n=3202 phased panel (GRCh38): resources/scripts/process_1kGP_3202_panel.sh.
#   GENCODE v44: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
#   Prebuilt genome-wide hg38 window BED + blacklist (resources/data/) are used directly.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../../.." && pwd)"
OUT="${1:-${REPO}/.test-run/COLO829_full/reference}"
mkdir -p "${OUT}"

UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
GENCODE_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"

CHROMS=$(seq 1 22)

echo "[1/3] reference FASTA (hg38 chr1-chr22)"
if [ ! -f "${OUT}/hg38.fa.fai" ]; then
  : > "${OUT}/hg38.fa"
  for i in ${CHROMS}; do
    curl -fsSL "${UCSC}/chr${i}.fa.gz" | gunzip -c >> "${OUT}/hg38.fa"
  done
  samtools faidx "${OUT}/hg38.fa"
fi

echo "[2/3] GENCODE v44 GTF (chr1-chr22)"
if [ ! -f "${OUT}/gencode.hg38.gtf.gz" ]; then
  curl -fsSL "${GENCODE_GTF}" | gunzip -c \
    | awk -F'\t' '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/' \
    | gzip -c > "${OUT}/gencode.hg38.gtf.gz"
fi

echo "[3/3] 1000G n=3202 het-SNP targets (chr1-chr22)"
if [ ! -f "${OUT}/target_positions/target.chr22.pos.gz.tbi" ]; then
  bash "${REPO}/resources/scripts/process_1kGP_3202_panel.sh" \
    --ref hg38 --chroms "${CHROMS}" "${OUT}"
fi

echo "done: staged COLO829 full-autosome reference under ${OUT}"
