#!/usr/bin/env bash
set -euo pipefail

# Prepare 1kGP n=3,202 phased panel for use with the genotyping pipeline.
# Produces per-chromosome SNP panel VCFs, target position files, and phasing
# panel BCFs.
#
# Supports two reference versions:
#   hg38     — downloads per-chromosome VCFs from EBI FTP
#   chm13v2  — downloads (or uses existing) whole-genome BCF from HPRC S3
#
# Requires: bcftools, bgzip, tabix, wget (for download)

HG38_FTP="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
CHM13V2_S3="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/Phased_SHAPEIT5_v1.1"
CHM13V2_BCF="1KGP.CHM13v2.0.whole_genome.recalibrated.snp_indel.pass.phased.native_maps.biallelic.3202.bcf.gz"

usage() {
  cat >&2 <<EOF
Usage: $0 --ref <hg38|chm13v2> <output_dir> [--bcf /path/to/existing.bcf.gz] [--chroms "8 10"]

  --ref      Reference version (required): hg38 or chm13v2
  output_dir Directory to write outputs
  --bcf      Path to pre-downloaded whole-genome BCF (chm13v2 only; skips download)
  --chroms   Space-separated chromosomes to process (default: "1..22 X"); use to
             stage only a subset, e.g. --chroms "8 10" for a chr8/chr10 test case
EOF
  exit 1
}

REFVERS=""
OUT=""
BCF_INPUT=""
CHROMS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ref) REFVERS="$2"; shift 2 ;;
    --bcf) BCF_INPUT="$2"; shift 2 ;;
    --chroms) CHROMS="$2"; shift 2 ;;
    -h|--help) usage ;;
    *)
      if [[ -z "${OUT}" ]]; then
        OUT="$1"; shift
      else
        usage
      fi
      ;;
  esac
done

[[ -z "${REFVERS}" || -z "${OUT}" ]] && usage
[[ "${REFVERS}" != "hg38" && "${REFVERS}" != "chm13v2" ]] && { echo "ERROR: --ref must be hg38 or chm13v2" >&2; exit 1; }
[[ -z "${CHROMS}" ]] && CHROMS="$(seq 1 22) X"

mkdir -p "${OUT}"/{snps,target_positions,phasing_panel}

##################################################
# Download / locate source data
##################################################

if [[ "${REFVERS}" == "hg38" ]]; then
  [[ -n "${BCF_INPUT}" ]] && { echo "WARNING: --bcf ignored for hg38 (per-chromosome download)" >&2; }
  mkdir -p "${OUT}/raw"
fi

if [[ "${REFVERS}" == "chm13v2" ]]; then
  if [[ -n "${BCF_INPUT}" ]]; then
    [[ ! -f "${BCF_INPUT}" ]] && { echo "ERROR: BCF not found: ${BCF_INPUT}" >&2; exit 1; }
    echo "using existing BCF: ${BCF_INPUT}"
    INPUT_BCF="${BCF_INPUT}"
  else
    mkdir -p "${OUT}/raw"
    INPUT_BCF="${OUT}/raw/${CHM13V2_BCF}"
    if [ ! -f "${INPUT_BCF}" ]; then
      date; echo "downloading ${CHM13V2_BCF}"
      wget -c -P "${OUT}/raw" "${CHM13V2_S3}/${CHM13V2_BCF}"
      wget -c -P "${OUT}/raw" "${CHM13V2_S3}/${CHM13V2_BCF}.csi"
    fi
  fi
fi

##################################################
# Per-chromosome processing
##################################################

: > "${OUT}/snp_vcfs.lst"
for chr in ${CHROMS}; do
  date
  echo "chr${chr}"

  if [[ "${REFVERS}" == "hg38" ]]; then
    # hg38: download per-chromosome VCF from EBI FTP
    if [[ "${chr}" == "X" ]]; then
      vcfgz="1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    else
      vcfgz="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    fi
    if [ ! -f "${OUT}/raw/${vcfgz}" ]; then
      wget -c -P "${OUT}/raw" "${HG38_FTP}/${vcfgz}"
      wget -c -P "${OUT}/raw" "${HG38_FTP}/${vcfgz}.tbi"
    fi
    SRC_VCF="${OUT}/raw/${vcfgz}"
  else
    # chm13v2: source is the whole-genome BCF
    SRC_VCF="${INPUT_BCF}"
  fi

  # SNP sites (no genotypes)
  if [ ! -f "${OUT}/snps/chr${chr}.vcf.gz" ]; then
    echo "extract SNP sites"
    if [[ "${REFVERS}" == "hg38" ]]; then
      bcftools view -v snps -G -Oz -o "${OUT}/snps/chr${chr}.vcf.gz" "${SRC_VCF}"
    else
      bcftools view -v snps -G -r "chr${chr}" -Oz -o "${OUT}/snps/chr${chr}.vcf.gz" "${SRC_VCF}"
    fi
    tabix -p vcf "${OUT}/snps/chr${chr}.vcf.gz"
  fi
  echo "${OUT}/snps/chr${chr}.vcf.gz" >> "${OUT}/snp_vcfs.lst"

  # Target positions
  if [ ! -f "${OUT}/target_positions/target.chr${chr}.pos.gz" ]; then
    echo "extract target positions"
    bcftools query -f '%CHROM\t%POS\n' \
      "${OUT}/snps/chr${chr}.vcf.gz" | bgzip -c > "${OUT}/target_positions/target.chr${chr}.pos.gz"
    tabix -s1 -b2 -e2 "${OUT}/target_positions/target.chr${chr}.pos.gz"
  fi

  # Phasing panel BCF
  if [ ! -f "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" ]; then
    echo "extract phasing panel BCF"
    if [[ "${REFVERS}" == "hg38" ]]; then
      bcftools view -Ob -o "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" "${SRC_VCF}"
    else
      bcftools view -r "chr${chr}" -Ob -o "${OUT}/phasing_panel/chr${chr}.genotypes.bcf" "${SRC_VCF}"
    fi
    bcftools index --csi "${OUT}/phasing_panel/chr${chr}.genotypes.bcf"
  fi
done

# Merged SNP panel
bcftools concat -f "${OUT}/snp_vcfs.lst" -Oz -o "${OUT}/snps.vcf.gz"
tabix -p vcf "${OUT}/snps.vcf.gz"

##################################################
# Cleanup intermediates
##################################################
echo "=== Cleaning up intermediates ==="
rm -rf "${OUT}/raw" "${OUT}/snps" "${OUT}/snp_vcfs.lst" "${OUT}/snp_files.lst"

date
echo "=====summary====="
echo "reference: ${REFVERS}"
echo "target positions: ${OUT}/target_positions/"
echo "phasing panel: ${OUT}/phasing_panel/"
echo "snp panel: ${OUT}/snps.vcf.gz"
exit 0
