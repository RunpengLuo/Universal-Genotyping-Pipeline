#!/usr/bin/env bash
# Stage hg38 chr6+chr16 reference + 1000G panel + Eagle map for the HCC1395 test case.
#
# Runpeng Luo (2026-07-26)
#
# Dependencies (on PATH): curl, wget, samtools, bcftools, bgzip, tabix (+ awk/gzip/tar).
#   A minimal env suffices, e.g.
#   `mamba create -n genotyping-prep -c bioconda -c conda-forge bcftools samtools htslib wget`.
#
# Usage:
#   bash prepare_refs.sh [out_dir]
#           out_dir  # default: <repo>/.test-run/HCC1395/reference (gitignored)
# Outputs (match tests/data/HCC1395/config.yaml keys):
#   hg38.fa (+ .fai)                           -> reference (chr6 + chr16)
#   gencode.hg38.gtf.gz                         -> gtf_file (chr6 + chr16)
#   target_positions/target.chr{6,16}.pos.gz   -> snp_targets (1000G n=3202)
#   phasing_panel/chr{6,16}.genotypes.bcf       -> phasing_panel (Eagle reference)
#   genetic_map_hg38_withX.txt.gz               -> gmap_path (Eagle2)
# Notes/References:
#   Reference: UCSC hg38 per-chromosome FASTA (GRCh38 primary; matches the chr-prefixed
#     GRCh38.d1.vd1 Illumina BAMs). Statistical phasing (Eagle2) against the 1000G panel.
#   1000G n=3202 phased panel (GRCh38): resources/scripts/process_1kGP_3202_panel.sh.
#   Eagle2: https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
#   GENCODE v44: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
#   Prebuilt genome-wide hg38 window BED + blacklist (resources/data/) are used directly.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../../.." && pwd)"
OUT="${1:-${REPO}/.test-run/HCC1395/reference}"
TMP="${OUT}/_tmp"
trap 'rm -rf "${TMP}"' EXIT
mkdir -p "${OUT}" "${TMP}"

UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
GENCODE_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
EAGLE_TAR="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz"

echo "[1/4] reference FASTA (hg38 chr6 + chr16)"
if [ ! -f "${OUT}/hg38.fa.fai" ]; then
  : > "${OUT}/hg38.fa"
  for c in chr6 chr16; do
    curl -fsSL "${UCSC}/${c}.fa.gz" | gunzip -c >> "${OUT}/hg38.fa"
  done
  samtools faidx "${OUT}/hg38.fa"
fi

echo "[2/4] GENCODE v44 GTF (chr6 + chr16)"
if [ ! -f "${OUT}/gencode.hg38.gtf.gz" ]; then
  curl -fsSL "${GENCODE_GTF}" | gunzip -c \
    | awk -F'\t' '$1=="chr6" || $1=="chr16"' \
    | gzip -c > "${OUT}/gencode.hg38.gtf.gz"
fi

echo "[3/4] 1000G n=3202 targets + phasing panel (chr6, chr16)"
if [ ! -f "${OUT}/phasing_panel/chr16.genotypes.bcf.csi" ]; then
  bash "${REPO}/resources/scripts/process_1kGP_3202_panel.sh" \
    --ref hg38 --chroms "6 16" "${OUT}"
fi

echo "[4/4] Eagle hg38 genetic map"
if [ ! -f "${OUT}/genetic_map_hg38_withX.txt.gz" ]; then
  curl -fsSL "${EAGLE_TAR}" -o "${TMP}/eagle.tar.gz"
  tar -xzf "${TMP}/eagle.tar.gz" -C "${TMP}"
  cp "${TMP}"/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz "${OUT}/"
fi

rm -rf "${TMP}"
echo "done: staged HCC1395 reference under ${OUT}"
