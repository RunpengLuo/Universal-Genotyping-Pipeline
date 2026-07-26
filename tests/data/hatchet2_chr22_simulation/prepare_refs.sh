#!/usr/bin/env bash
# Stage hg19 chr22 reference + panel + tracks for hatchet2_chr22_simulation.
#
# Runpeng Luo (2026-07-25)
#
# Dependencies (on PATH): curl, bcftools, samtools, bgzip, tabix (+ coreutils awk/paste
#   /gzip/tar). A minimal env suffices, e.g.
#   `mamba create -n genotyping-prep -c bioconda -c conda-forge bcftools samtools htslib`.
#
# Usage:
#   bash prepare_refs.sh [out_dir]
#           out_dir   # where to stage files; default: <repo>/.test-run/<case>/reference
#                     # (gitignored runtime dir; config.yaml points here)
#
# Inputs (downloaded):
#   UCSC hg19 chr22 FASTA; HATCHet CI phasing panel (Zenodo 6709541, 1000GP Phase3
#   chr22 hap/legend/sample, GRCh37); Eagle v2.4.1 genetic map; GENCODE v19 GTF;
#   ENCODE hg19 blacklist v2.
# Outputs (match tests/.../config.yaml keys):
#   hg19.chr22.fa (+ .fai)                     -> reference
#   phasing_panel/chr22.genotypes.bcf (+ .csi) -> phasing_panel (biallelic SNPs, MAF>=MAF)
#   snp_targets/target.chr22.pos.gz (+ .tbi)   -> snp_targets (dir)
#   genetic_map_hg19_withX.txt.gz              -> gmap_path
#   gencode.v19.chr22.gtf.gz                   -> gtf_file
#   hg19-blacklist.v2.chr22.bed.gz             -> blacklist_bed
#   windows.1kbp.hg19.chr22.bed.gz             -> window_bed (pre-built, chr22 subset)
# Notes/References:
#   Phasing panel: https://zenodo.org/records/6709541 (HATCHet CI, 1000GP Phase3 chr22
#     in IMPUTE2 hap/legend/sample; converted to a chr-prefixed BCF for Eagle here).
#   UCSC hg19: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/
#   Eagle2: https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
#   GENCODE v19: https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/
#   ENCODE blacklist: https://github.com/Boyle-Lab/Blacklist/tree/master/lists
#   The pre-built window BED (resources/data/windows.1kbp.hg19.bed.gz) already carries
#   GC/MAP/REPLI, so window_bed skips build_window_bed + Repli-seq (no mappability track
#   needed). All tracks use chr-prefixed contigs (chr22), matching the hg19 BAMs.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "${HERE}/../../.." && pwd)"
OUT="${1:-${REPO}/.test-run/hatchet2_chr22_simulation/reference}"
TMP="${OUT}/_tmp"
trap 'rm -rf "${TMP}"' EXIT
mkdir -p "${OUT}/snp_targets" "${OUT}/phasing_panel" "${TMP}"

MAF=0.01  # minor-allele-freq cut: keep all 2504 samples, common SNPs only (bounds eagle memory)

UCSC_FA="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz"
ZENODO_PANEL="https://zenodo.org/records/6709541/files/1000GP_Phase3.tgz?download=1"
EAGLE_TAR="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz"
GENCODE_GTF="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
BLACKLIST="https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg19-blacklist.v2.bed.gz"

echo "[1/6] reference FASTA (hg19 chr22)"
if [ ! -f "${OUT}/hg19.chr22.fa.fai" ]; then
  curl -fsSL "${UCSC_FA}" | gunzip -c > "${OUT}/hg19.chr22.fa"
  samtools faidx "${OUT}/hg19.chr22.fa"
fi

echo "[2/6] phasing panel + SHAPEIT map (Zenodo 6709541 hap/legend -> BCF, MAF>=${MAF})"
PDIR="${TMP}/1000GP_Phase3"
need_panel=0; [ -f "${OUT}/phasing_panel/chr22.genotypes.bcf.csi" ] || need_panel=1
need_map=0;   [ -f "${OUT}/maps/chr22.b37.gmap.gz" ] || need_map=1
if [ "${need_panel}" = 1 ] || [ "${need_map}" = 1 ]; then
  curl -fsSL "${ZENODO_PANEL}" -o "${TMP}/1000GP_Phase3.tgz"
  tar -xzf "${TMP}/1000GP_Phase3.tgz" -C "${TMP}"
fi
if [ "${need_panel}" = 1 ]; then
  # IMPUTE2 .samples: 2-line header (ID_1 ID_2 missing / 0 0 0) then "id id 0"
  {
    echo "ID_1 ID_2 missing"
    echo "0 0 0"
    tail -n +2 "${PDIR}/1000GP_Phase3.sample" | awk '{print $1, $1, 0}'
  } > "${TMP}/panel.samples"
  # ShapeIt .haps: leading cols CHROM, CHROM:POS_REF_ALT, POS, a0, a1 (row-aligned with
  # the hap matrix); col2 is where bcftools reads CHROM, so chr-prefix it here.
  paste -d' ' \
    <(gunzip -c "${PDIR}/1000GP_Phase3_chr22.legend.gz" | tail -n +2 \
        | awk '{print "chr22", "chr22:"$2"_"$3"_"$4, $2, $3, $4}') \
    <(gunzip -c "${PDIR}/1000GP_Phase3_chr22.hap.gz") \
    | gzip -c > "${TMP}/panel.haps.gz"
  bcftools convert --hapsample2vcf "${TMP}/panel.haps.gz","${TMP}/panel.samples" -Ou \
    | bcftools view -v snps -m2 -M2 -q "${MAF}:minor" \
        -Ob -o "${OUT}/phasing_panel/chr22.genotypes.bcf"
  bcftools index --csi "${OUT}/phasing_panel/chr22.genotypes.bcf"
fi
if [ "${need_map}" = 1 ]; then
  # SHAPEIT5 --map: "pos chr cM"; chr column matches the chr-prefixed --region (chr22).
  # Source is the 1000GP combined map (position COMBINED_rate Genetic_Map(cM)).
  mkdir -p "${OUT}/maps"
  {
    printf "pos\tchr\tcM\n"
    tail -n +2 "${PDIR}/genetic_map_chr22_combined_b37.txt" \
      | awk '{print $1"\tchr22\t"$3}'
  } | gzip -c > "${OUT}/maps/chr22.b37.gmap.gz"
fi
if [ ! -f "${OUT}/snp_targets/target.chr22.pos.gz.tbi" ]; then
  bcftools query -f '%CHROM\t%POS\n' "${OUT}/phasing_panel/chr22.genotypes.bcf" \
    | bgzip -c > "${OUT}/snp_targets/target.chr22.pos.gz"
  tabix -s1 -b2 -e2 "${OUT}/snp_targets/target.chr22.pos.gz"
fi

echo "[3/6] Eagle hg19 genetic map"
if [ ! -f "${OUT}/genetic_map_hg19_withX.txt.gz" ]; then
  curl -fsSL "${EAGLE_TAR}" -o "${TMP}/eagle.tar.gz"
  tar -xzf "${TMP}/eagle.tar.gz" -C "${TMP}"
  cp "${TMP}"/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz "${OUT}/"
fi

echo "[4/6] GENCODE v19 GTF (chr22)"
if [ ! -f "${OUT}/gencode.v19.chr22.gtf.gz" ]; then
  curl -fsSL "${GENCODE_GTF}" | gunzip -c | awk -F'\t' '$1=="chr22"' \
    | gzip -c > "${OUT}/gencode.v19.chr22.gtf.gz"
fi

echo "[5/6] ENCODE hg19 blacklist v2 (chr22)"
if [ ! -f "${OUT}/hg19-blacklist.v2.chr22.bed.gz" ]; then
  curl -fsSL "${BLACKLIST}" | gunzip -c | awk -F'\t' '$1=="chr22"' \
    | bgzip -c > "${OUT}/hg19-blacklist.v2.chr22.bed.gz"
fi

echo "[6/6] pre-built hg19 window BED, subset to chr22"
if [ ! -f "${OUT}/windows.1kbp.hg19.chr22.bed.gz" ]; then
  gzip -dc "${REPO}/resources/data/windows.1kbp.hg19.bed.gz" \
    | awk 'NR==1 || $1=="chr22"' \
    | gzip -c > "${OUT}/windows.1kbp.hg19.chr22.bed.gz"
fi

rm -rf "${TMP}"
echo "done: staged reference under ${OUT}"
