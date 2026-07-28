#!/usr/bin/env bash
set -euo pipefail

# Build a mouse SNP panel + phasing panel from the UCSC MGP v5 strain VCF.
#
# Produces three artifacts compatible with the universal-genotyping pipeline:
#   1. snp_panel        — sites-only biallelic SNP VCF (no GT columns).
#   2. phasing_panel/   — per-chromosome multi-strain BCFs with phased GTs.
#                         Inbred strain calls are phased trivially: 0/0->0|0,
#                         1/1->1|1, strain-het -> ./. (assumed noise in inbred).
#   3. target_positions — per-chromosome SNP position files for bcftools mpileup.
#
# Strains: by default all strains in the MGP v5 VCF (36) are used. Pass --strains
# to restrict to a subset (e.g. a B6 x 129 cross): both panels are subset to those
# strains and only sites polymorphic among them are kept, so every target is
# informative for that cross. A subset also shrinks the phasing reference, so keep
# the default when you want a diverse haplotype panel for eagle/shapeit.
#
# Chromosome naming: chr1..chr19,chrX (matches 10x mm10-2020-A reference).
#
# Output (under <out-dir>):
#   mgpV5.biallelic_snps.vcf.gz[.tbi]
#   phasing_panel/chr{1..19,X}.genotypes.bcf[.csi]
#   target_positions/target.chr{1..19,X}.pos.gz[.tbi]
#
# Usage:
#   bash build_mouse_mgp_panel.sh \
#       --out-dir /path/to/mouse_panel \
#       [--strains C57BL_6NJ,129S1_SvImJ] \
#       [--build-snp-targets /path/to/build_snp_targets.sh]
#
# Requires: bcftools, bgzip, tabix, wget, awk

MGP_URL="https://hgdownload.soe.ucsc.edu/gbdb/mm10/mouseStrains/mgpV5MergedSNPsAlldbSNP142.vcf.gz"
MGP_NAME="mgpV5MergedSNPsAlldbSNP142.vcf.gz"
MOUSE_CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X"

usage() {
  cat >&2 <<EOF
Usage: $0 --out-dir DIR [options]

  --out-dir DIR             Output root (required)
  --strains LIST            Comma-separated strain subset (default: all strains in the VCF).
                            Restricts both panels to these strains and to sites polymorphic
                            among them, e.g. C57BL_6NJ,129S1_SvImJ for a B6 x 129 cross.
  --build-snp-targets PATH  Path to build_snp_targets.sh. [default: sibling in this scripts/ dir]
EOF
  exit 1
}

OUT=""
STRAINS=""
BUILD_SNP_TARGETS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT="$2"; shift 2 ;;
    --strains) STRAINS="$2"; shift 2 ;;
    --build-snp-targets) BUILD_SNP_TARGETS="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown arg: $1" >&2; usage ;;
  esac
done

[[ -z "${OUT}" ]] && usage
mkdir -p "${OUT}/raw" "${OUT}/phasing_panel"

# strain subset: restrict samples (bcftools -s updates AC/AN) and keep only sites
# polymorphic within them (-c 1:minor). Empty by default -> all strains, no site filter.
STRAIN_ARG=""
POLY_ARG=""
if [[ -n "${STRAINS}" ]]; then
  STRAIN_ARG="-s ${STRAINS}"
  POLY_ARG="-c 1:minor"
fi

##################################################
# 1. Download MGP VCF
##################################################
MGP_VCF="${OUT}/raw/${MGP_NAME}"
if [[ ! -f "${MGP_VCF}" ]]; then
  date; echo "downloading ${MGP_NAME}"
  wget -c -P "${OUT}/raw" "${MGP_URL}"
  wget -c -P "${OUT}/raw" "${MGP_URL}.tbi"
fi

##################################################
# 2. Sanity: contigs + strain list
##################################################
echo "=== MGP VCF contigs (expect chr1..chr19,chrX) ==="
bcftools view -h "${MGP_VCF}" | grep ^##contig | head -25
echo
N_STRAINS=$(bcftools query -l "${MGP_VCF}" | wc -l | tr -d ' ')
echo "=== strains in MGP VCF: ${N_STRAINS} (using: ${STRAINS:-all}) ==="
bcftools query -l "${MGP_VCF}" | sed 's/^/  /'
echo

REGIONS=$(echo "${MOUSE_CHROMS}" | sed 's/ /,chr/g; s/^/chr/')
FILTER='FILTER="PASS" && TYPE="snp" && N_ALT=1'

##################################################
# 3. snp_panel: sites-only biallelic SNPs
##################################################
PANEL_VCF="${OUT}/mgpV5.biallelic_snps.vcf.gz"
if [[ ! -f "${PANEL_VCF}" ]]; then
  echo "[snp_panel] extracting biallelic SNP sites -> ${PANEL_VCF}"
  bcftools view "${MGP_VCF}" -r "${REGIONS}" -i "${FILTER}" ${STRAIN_ARG} ${POLY_ARG} \
      -G -Oz -o "${PANEL_VCF}"
  tabix -f -p vcf "${PANEL_VCF}"
fi

##################################################
# 4. phasing_panel: per-chrom multi-strain BCFs with phased GTs
##################################################
STATS_TSV="${OUT}/phasing_panel_stats.tsv"
printf "chr\tn_sites\tn_strain_het_masked\n" > "${STATS_TSV}"
for chr in ${MOUSE_CHROMS}; do
  bcf_out="${OUT}/phasing_panel/chr${chr}.genotypes.bcf"
  if [[ -f "${bcf_out}" && -f "${bcf_out}.csi" ]]; then
    echo "[phasing_panel] chr${chr} already exists, skip"
    continue
  fi
  date; echo "[phasing_panel] chr${chr}"
  bcftools view "${MGP_VCF}" -r "chr${chr}" -i "${FILTER}" ${STRAIN_ARG} ${POLY_ARG} \
    | bcftools annotate -x INFO,^FORMAT/GT \
    | awk -v chr="${chr}" -v stats="${STATS_TSV}" '
        BEGIN{OFS="\t"; n=0; mask=0}
        /^#/ {print; next}
        {
          n++; row_had_het=0
          for (i=10; i<=NF; i++) {
            g = $i
            if      (g == "0/0") $i = "0|0"
            else if (g == "1/1") $i = "1|1"
            else if (g == "0/1" || g == "1/0") {$i = "./."; row_had_het=1}
          }
          if (row_had_het) mask++
          print
        }
        END{ printf "chr%s\t%d\t%d\n", chr, n, mask >> stats }
      ' \
    | bcftools view -Ob -o "${bcf_out}"
  bcftools index --csi -f "${bcf_out}"
done

##################################################
# 5. target_positions/ for bcftools mpileup -T
##################################################
if [[ -z "${BUILD_SNP_TARGETS}" ]]; then
  SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
  CANDIDATE="${SCRIPT_DIR}/build_snp_targets.sh"
  if [[ -f "${CANDIDATE}" ]]; then
    BUILD_SNP_TARGETS="${CANDIDATE}"
  else
    echo "ERROR: build_snp_targets.sh not found; pass --build-snp-targets PATH" >&2
    exit 1
  fi
fi
echo "[snp_targets] running build_snp_targets.sh"
bash "${BUILD_SNP_TARGETS}" --chroms "${MOUSE_CHROMS}" \
    "${PANEL_VCF}" "${OUT}/target_positions"

##################################################
# 6. Stats: multi-allelic SNPs dropped + final counts
##################################################
echo
echo "=== panel stats ==="
N_BIALLELIC=$(bcftools view -H "${PANEL_VCF}" | wc -l | tr -d ' ')
date; echo "scanning source for total SNP count on target chroms (PASS, any N_ALT)..."
N_ALL_SNPS=$(bcftools view -H -r "${REGIONS}" -i 'FILTER="PASS" && TYPE="snp"' ${STRAIN_ARG} ${POLY_ARG} "${MGP_VCF}" | wc -l | tr -d ' ')
N_MULTI=$((N_ALL_SNPS - N_BIALLELIC))
PCT_MULTI=$(awk -v m="${N_MULTI}" -v t="${N_ALL_SNPS}" 'BEGIN{printf "%.2f", (t>0 ? 100*m/t : 0)}')

printf "strains:                   %s\n" "${STRAINS:-all (${N_STRAINS})}"
printf "all PASS SNPs (any N_ALT): %s\n" "${N_ALL_SNPS}"
printf "biallelic kept:            %s\n" "${N_BIALLELIC}"
printf "multi-allelic dropped:     %s (%s%%)\n" "${N_MULTI}" "${PCT_MULTI}"
echo
echo "per-chrom phasing_panel (sites / rows with any strain-het masked):"
column -t -s $'\t' "${STATS_TSV}" | sed 's/^/  /'

##################################################
# 7. Cleanup raw download
##################################################
rm -rf "${OUT}/raw"

date
echo "=====summary====="
echo "snp_panel:      ${PANEL_VCF}"
echo "phasing_panel:  ${OUT}/phasing_panel/   ($(ls "${OUT}/phasing_panel"/*.bcf 2>/dev/null | wc -l | tr -d ' ') chrom BCFs)"
echo "snp_targets:    ${OUT}/target_positions/"
echo "stats TSV:      ${STATS_TSV}"
exit 0
