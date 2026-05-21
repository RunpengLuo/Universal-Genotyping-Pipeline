#!/usr/bin/env bash
set -euo pipefail

# Build per-chromosome SNP target position files from any SNP panel VCF.
# Output: target.chr<C>.pos.gz + .tbi in the specified output directory.
#
# Usage:
#   bash build_snp_targets.sh [--chroms "1 2 ... X"] <input.vcf.gz> <output_dir>
#
# --chroms : space-separated list of chromosome short names (without "chr" prefix).
#            Default: "1 2 3 ... 22 X" (human).
#            For mouse pass: --chroms "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X"
#
# Requires: bcftools, bgzip, tabix

usage() {
  echo "Usage: $0 [--chroms \"1 2 ... X\"] <input.vcf.gz> <output_dir>" >&2
  exit 1
}

CHROMS_STR=""
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --chroms) CHROMS_STR="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) POSITIONAL+=("$1"); shift ;;
  esac
done

[[ ${#POSITIONAL[@]} -eq 2 ]] || usage

INPUT_VCF=${POSITIONAL[0]}
OUT=${POSITIONAL[1]}

if [[ -z "${CHROMS_STR}" ]]; then
  CHROMS_STR="$(seq 1 22) X"
fi
read -r -a CHROMS <<< "${CHROMS_STR}"

if [[ ! -f "${INPUT_VCF}" ]]; then
  echo "Error: input VCF not found: ${INPUT_VCF}" >&2
  exit 1
fi

mkdir -p "${OUT}"

for chr in "${CHROMS[@]}"; do
  outfile="${OUT}/target.chr${chr}.pos.gz"
  if [[ -f "${outfile}" && -f "${outfile}.tbi" ]]; then
    echo "skip chr${chr} (already exists)"
    continue
  fi
  echo "chr${chr}"
  bcftools query -r "chr${chr}" -f '%CHROM\t%POS\n' "${INPUT_VCF}" \
    | bgzip -c > "${outfile}"
  tabix -s1 -b2 -e2 "${outfile}"
done

echo "Done. Output: ${OUT}/target.chr{${CHROMS[*]}}.pos.gz"
