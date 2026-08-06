#!/usr/bin/env bash
# Run cellranger-atac count on one or more scATAC samples, logging each run.
#
# Runpeng Luo
# Last update: 2026-08-06
#
# Dependencies:
#   cellranger-atac (set CELLRANGER); a Cell Ranger ARC / ATAC reference (set REFERENCE).
#
# Usage:
#   Edit the CONFIG block and the SAMPLES list below, then:
#     bash run_cellranger_atac.sh
#   Or override any config value from the environment:
#     CELLRANGER=/path/cellranger-atac REFERENCE=/path/ref WORKDIR=/path/out \
#       FASTQ_ROOT=/path/fastq LOCALCORES=32 LOCALMEM=256 bash run_cellranger_atac.sh
#
# Inputs
#   REFERENCE:  cellranger-atac/-arc reference directory
#   FASTQ_ROOT: parent of the per-sample fastq directories named in SAMPLES
# Outputs:
#   WORKDIR/<id>/outs/ per sample; WORKDIR/logs/<id>_<timestamp>.log per run
# Notes/References:
#   cellranger-atac count: https://www.10xgenomics.com/support/software/cell-ranger-atac
#   --fastqs is a directory holding <sample>_S*_L*_R{1,2,3}_*.fastq.gz for one sample.

set -uo pipefail

# ---- CONFIG (override via environment) -------------------------------------
CELLRANGER="${CELLRANGER:?set CELLRANGER to the cellranger-atac binary}"
REFERENCE="${REFERENCE:?set REFERENCE to a cellranger-atac/-arc reference dir}"
WORKDIR="${WORKDIR:-$(pwd)}"
FASTQ_ROOT="${FASTQ_ROOT:?set FASTQ_ROOT to the parent of the fastq dirs}"
LOCALCORES="${LOCALCORES:-16}"
LOCALMEM="${LOCALMEM:-128}"

# ---- SAMPLES: "<output_id> <fastq_subdir>", one per line --------------------
# <fastq_subdir> is relative to FASTQ_ROOT and holds one sample's fastqs.
SAMPLES=(
    "sample_name1  fastq_dir1"
    "sample_name2  fastq_dir2"
)

LOG_DIR="${WORKDIR}/logs"

run_cr_atac() {
    local id=$1
    local fastqs=$2
    local timestamp logfile status

    timestamp=$(date '+%Y%m%d_%H%M%S')
    logfile="${LOG_DIR}/${id}_${timestamp}.log"

    {
        echo "============================================================"
        echo "Sample:     ${id}"
        echo "FASTQs:     ${fastqs}"
        echo "Reference:  ${REFERENCE}"
        echo "Started:    $(date)"
        echo "Log:        ${logfile}"
        echo "============================================================"
    } | tee "${logfile}"

    "${CELLRANGER}" count \
        --id="${id}" \
        --reference="${REFERENCE}" \
        --fastqs="${fastqs}" \
        --localcores="${LOCALCORES}" \
        --localmem="${LOCALMEM}" \
        2>&1 | tee -a "${logfile}"

    status=${PIPESTATUS[0]}

    if (( status == 0 )); then
        echo "SUCCESS: ${id} completed at $(date)" | tee -a "${logfile}"
    else
        echo "ERROR: ${id} failed with exit status ${status} at $(date)" \
            | tee -a "${logfile}" >&2
    fi

    return 0
}

main() {
    cd "${WORKDIR}"
    mkdir -p "${LOG_DIR}"

    local failed=0
    for entry in "${SAMPLES[@]}"; do
        read -r id subdir <<<"${entry}"
        run_cr_atac "${id}" "${FASTQ_ROOT}/${subdir}" || failed=1
    done
    return "${failed}"
}

main "$@"
