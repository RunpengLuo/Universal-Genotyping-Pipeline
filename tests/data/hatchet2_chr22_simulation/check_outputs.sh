#!/usr/bin/env bash
# Assert a tests/data/<case> produced the expected bulk outputs. A step of the CI job.
#
# Runpeng Luo (2026-07-24)
#
# Usage:
#   [RUN=<out_dir>] bash check_outputs.sh
#           RUN   # output root (default: <repo>/.test-run/<case>)
# Outputs checked (per min_snp_reads value): bb.tsv.gz, bb.{Tallele,Aallele,Ballele,
#   depth,rdr}.npz, sample_ids.tsv under RUN/bb/MSR{msr}/bulk/. Exits nonzero on any
#   missing/empty file or an empty bin table. CASE is inferred from this script's
#   location, so the file is identical across cases.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
RUN="${RUN:-$REPO/.test-run/$(basename "$HERE")}"
# Auto-detect the MSR{msr} dirs the run produced (matches whatever min_snp_reads was
# configured); override with MSRS="500 1000" to pin specific values.
MSRS="${MSRS:-$(ls -d "$RUN"/bb/MSR*/ 2>/dev/null | sed 's#.*/MSR##; s#/##' | tr '\n' ' ')}"
[ -n "${MSRS// /}" ] || { echo "no bb/MSR*/ dirs under $RUN"; exit 1; }

fail=0
for msr in $MSRS; do
  d="$RUN/bb/MSR${msr}/bulk"
  for f in bb.tsv.gz bb.Tallele.npz bb.Aallele.npz bb.Ballele.npz \
           bb.depth.npz bb.rdr.npz sample_ids.tsv; do
    if [ ! -s "$d/$f" ]; then echo "MISSING/EMPTY: $d/$f"; fail=1; fi
  done
done
[ "$fail" -eq 0 ] || { echo "output check FAILED"; exit 1; }

for msr in $MSRS; do
  bb="$RUN/bb/MSR${msr}/bulk/bb.tsv.gz"
  echo "=== bb.tsv.gz (MSR${msr}) head ==="
  head -n 10 <(zcat "$bb")   # <(...) avoids zcat|head SIGPIPE under pipefail
  nbins=$(( $(zcat "$bb" | wc -l) - 1 ))
  echo "bins (MSR${msr}): $nbins"
  [ "$nbins" -gt 0 ] || { echo "no bins produced for MSR${msr}"; exit 1; }
done
echo "output check PASSED"
