#!/usr/bin/env bash
# Stage hg38 reference + 1000G panel + Eagle map for the HCC1395 test case.
#
# Runpeng Luo (2026-07-26)
#
# Dependencies (on PATH): curl, wget, samtools, bcftools, bgzip, tabix (+ awk/gzip/tar).
#   A minimal env suffices, e.g.
#   `mamba create -n genotyping-prep -c bioconda -c conda-forge bcftools samtools htslib wget`.
#
# Usage:
#   [CHROMS="1 2 ... 22"] bash prepare_refs.sh [out_dir]
#           CHROMS   # space-separated chromosome numbers to stage (default: "16").
#                    #   chr16-only for the default CI config; use "$(seq 1 22)" for
#                    #   the whole-genome run (run_pipeline.whole_genome.sh).
#           out_dir  # default: <repo>/.test-run/HCC1395/reference (gitignored)
# Outputs (match tests/data/HCC1395/config.yaml keys), for each c in CHROMS:
#   hg38.fa (+ .fai)                           -> reference
#   gencode.hg38.gtf.gz                         -> gtf_file
#   target_positions/target.chr${c}.pos.gz      -> snp_targets (1000G n=3202)
#   phasing_panel/chr${c}.genotypes.bcf         -> phasing_panel (Eagle reference)
#   genetic_map_hg38_withX.txt.gz               -> gmap_path (Eagle2, genome-wide)
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
CHROMS="${CHROMS:-16}"
LAST="${CHROMS##* }"
TMP="${OUT}/_tmp"
trap 'rm -rf "${TMP}"' EXIT
mkdir -p "${OUT}" "${TMP}"

# Restage the chrom-scoped files if CHROMS changed since the last run (the gmap is
# genome-wide, so it is not scoped). The per-step guards below then rebuild them.
STAMP="${OUT}/.staged_chroms"
if [ -f "${STAMP}" ] && [ "$(cat "${STAMP}")" != "${CHROMS}" ]; then
  echo "chrom scope changed ($(cat "${STAMP}") -> ${CHROMS}); restaging"
  rm -f "${OUT}/hg38.fa" "${OUT}/hg38.fa.fai" "${OUT}/gencode.hg38.gtf.gz"
  rm -rf "${OUT}/phasing_panel" "${OUT}/target_positions"
fi

UCSC="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes"
GENCODE_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
EAGLE_TAR="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz"

echo "[1/4] reference FASTA (hg38 chr: ${CHROMS})"
if [ ! -f "${OUT}/hg38.fa.fai" ]; then
  : > "${OUT}/hg38.fa"
  for c in ${CHROMS}; do
    curl -fsSL "${UCSC}/chr${c}.fa.gz" | gunzip -c >> "${OUT}/hg38.fa"
  done
  samtools faidx "${OUT}/hg38.fa"
fi

echo "[2/4] GENCODE v44 GTF (chr: ${CHROMS})"
if [ ! -f "${OUT}/gencode.hg38.gtf.gz" ]; then
  curl -fsSL "${GENCODE_GTF}" | gunzip -c \
    | awk -F'\t' -v cl="${CHROMS}" \
        'BEGIN{n=split(cl,a," ");for(i=1;i<=n;i++)keep["chr"a[i]]=1} ($1 in keep)' \
    | gzip -c > "${OUT}/gencode.hg38.gtf.gz"
fi

echo "[3/4] 1000G n=3202 targets + phasing panel (chr: ${CHROMS})"
if [ ! -f "${OUT}/phasing_panel/chr${LAST}.genotypes.bcf.csi" ]; then
  bash "${REPO}/resources/scripts/process_1kGP_3202_panel.sh" \
    --ref hg38 --chroms "${CHROMS}" "${OUT}"
fi

echo "[4/4] Eagle hg38 genetic map"
if [ ! -f "${OUT}/genetic_map_hg38_withX.txt.gz" ]; then
  curl -fsSL "${EAGLE_TAR}" -o "${TMP}/eagle.tar.gz"
  tar -xzf "${TMP}/eagle.tar.gz" -C "${TMP}"
  cp "${TMP}"/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz "${OUT}/"
fi

echo "${CHROMS}" > "${STAMP}"
rm -rf "${TMP}"
echo "done: staged HCC1395 reference (chr: ${CHROMS}) under ${OUT}"
