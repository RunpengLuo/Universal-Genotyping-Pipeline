#!/usr/bin/env bash
# Stage hg19 chr22 references for the dbSNP151 + SHAPEIT5 test variant.
#
# Runpeng Luo (2026-07-26)
#
# Runs the base prepare_refs.sh (shared refs + phasing panel + SHAPEIT map) into the
# dbSNP151 reference dir, then OVERRIDES snp_targets with dbSNP151 GRCh37 GATK chr22
# positions (chr-prefixed, so no renaming). Matches how the HATCHet2 demo data was
# simulated: germline variants seeded at dbSNP151 positions.
#
# Dependencies (on PATH): curl, bcftools, samtools, bgzip, tabix (+ awk/paste/gzip/tar).
#
# Usage:
#   bash prepare_refs.dbsnp151.sh [out_dir]
#           out_dir  # default: <repo>/.test-run-dbsnp151/<case>/reference
# Outputs (match config.dbsnp151.yaml keys):
#   snp_targets/target.chr22.pos.gz (+ .tbi)  -> snp_targets (dbSNP151 GATK chr22)
#   maps/chr22.b37.gmap.gz                     -> gmap_path (SHAPEIT: pos chr cM)
#   phasing_panel/, hg19.chr22.fa, blacklist, gtf, window bed  -> from prepare_refs.sh
# Notes/References:
#   dbSNP151 GRCh37 GATK build (chr-prefixed contigs), remote-sliced to chr22 via tabix:
#     https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz
#   HATCHet uses this exact file for hg19 + chr_notation (run.py snps_mapping).
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../../.." && pwd)"
OUT="${1:-${REPO}/.test-run-dbsnp151/hatchet2_chr22_simulation/reference}"

DBSNP_URL="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz"

echo "[dbsnp151] base refs + phasing panel + SHAPEIT map -> ${OUT}"
bash "${HERE}/prepare_refs.sh" "${OUT}"

echo "[dbsnp151] override snp_targets with dbSNP151 GATK chr22 (remote tabix slice)"
# .dbsnp151 marks the override done, so reruns skip the ~8M-position remote slice
# (base prepare_refs.sh builds a panel-derived snp_targets first; we replace it once).
if [ ! -f "${OUT}/snp_targets/.dbsnp151" ]; then
  rm -rf "${OUT}/snp_targets"
  mkdir -p "${OUT}/snp_targets"
  # 2-col CHROM<TAB>POS positions file (skip the VCF header), tabix-indexed for -T.
  tabix "${DBSNP_URL}" chr22 \
    | cut -f1,2 \
    | bgzip -c > "${OUT}/snp_targets/target.chr22.pos.gz"
  tabix -s1 -b2 -e2 "${OUT}/snp_targets/target.chr22.pos.gz"
  touch "${OUT}/snp_targets/.dbsnp151"
fi

NPOS=$(gzip -dc "${OUT}/snp_targets/target.chr22.pos.gz" | wc -l | tr -d ' ')
echo "done: dbSNP151 snp_targets = ${NPOS} chr22 positions under ${OUT}"
