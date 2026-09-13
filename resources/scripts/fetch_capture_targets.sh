#!/usr/bin/env bash
set -euo pipefail

# Fetch a hybrid-capture kit's target intervals from the UCSC exome-probeset track.
#
# The result is the `target_bed` config key: a plain BED3+ of the kit's targets, which
# splits the window grid into on- and off-target strata for bulkWES bias correction and
# RDR. Any BED3+ from the kit vendor works just as well; this only automates the UCSC
# copies, which are already on the pipeline's `chr` naming and lifted per build.
#
# Track: https://genome.ucsc.edu/cgi-bin/hgTrackUi?g=exomeProbesets
# Files: https://hgdownload.soe.ucsc.edu/gbdb/{hg19,hg38}/exomeProbesets/
#
# Common kit names (hg38), pass the stem without `.bb`:
#   xgen-exome-research-panel-targets-hg38        IDT xGen Exome Research Panel v1
#   xgen-exome-research-panel-v2-targets-hg38     IDT xGen Exome Research Panel v2
#   Twist_Exome_Target_hg38                       Twist Bioscience Exome
#   KAPA_HyperExome_hg38_primary_targets          Roche KAPA HyperExome
#   S07604514_Covered                             Agilent SureSelect Human All Exon V6
#   SeqCap_EZ_MedExome_hg38_capture_targets       Roche SeqCap EZ MedExome
#
# Dependencies:
#   curl, and bigBedToBed (conda: `bioconda::ucsc-bigbedtobed`)
#
# Usage:
#   bash fetch_capture_targets.sh <kit-stem> <build> <out-dir>
#           kit-stem   # file stem under the UCSC exomeProbesets directory, no .bb
#           build      # hg19 or hg38, must match the run's reference_version
#           out-dir    # directory for <kit-stem>.bed.gz
#
# Outputs:
#   <out-dir>/<kit-stem>.bed.gz: BED3 target intervals, chr-prefixed, 0-based half-open,
#     restricted to chr1-22/X/Y, coordinate-sorted, and overlapping intervals unioned

if [ "$#" -ne 3 ]; then
    sed -n '/^# Usage:/,/^$/p' "$0" >&2
    exit 1
fi

KIT="$1"
BUILD="$2"
OUT_DIR="$3"

case "$BUILD" in
    hg19 | hg38) ;;
    *)
        echo "ERROR: build must be hg19 or hg38, got '$BUILD'" >&2
        exit 1
        ;;
esac

command -v bigBedToBed >/dev/null || {
    echo "ERROR: bigBedToBed not on PATH (conda install -c bioconda ucsc-bigbedtobed)" >&2
    exit 1
}

URL="https://hgdownload.soe.ucsc.edu/gbdb/${BUILD}/exomeProbesets/${KIT}.bb"
mkdir -p "$OUT_DIR"
BB="${OUT_DIR}/${KIT}.bb"
RAW="${OUT_DIR}/${KIT}.raw.bed"
BED="${OUT_DIR}/${KIT}.bed.gz"

echo "[1/3] fetching ${URL}"
curl -fsS -o "$BB" "$URL"

echo "[2/3] converting to BED"
bigBedToBed "$BB" "$RAW"
rm -f "$BB"

# only columns 1-3 are ever read, and overlapping probes must count once, so normalize
# here: the pipeline unions internally either way, but a canonical file is smaller and
# reproducible. Alt contigs never appear in the window grid, so they are dropped.
echo "[3/3] normalizing (primary contigs, sorted, merged)"
N_ALT=$(awk '$1 !~ /^chr([0-9]+|X|Y)$/' "$RAW" | wc -l | tr -d ' ')
[ "$N_ALT" -gt 0 ] && echo "  dropping ${N_ALT} intervals on non-primary contigs"
awk '$1 ~ /^chr([0-9]+|X|Y)$/ {print $1 "\t" $2 "\t" $3}' "$RAW" \
    | LC_ALL=C sort -k1,1 -k2,2n \
    | awk 'BEGIN {OFS = "\t"}
           {
               if ($1 != c || $2 > e) { if (c != "") print c, s, e; c = $1; s = $2; e = $3 }
               else if ($3 > e) { e = $3 }
           }
           END { if (c != "") print c, s, e }' \
    | gzip -c > "$BED"
rm -f "$RAW"

echo "wrote ${BED}"
gzip -dc "$BED" \
    | awk '{n++; bp += $3 - $2} END {printf "  %d intervals, %.2f Mbp of target\n", n, bp / 1e6}'
