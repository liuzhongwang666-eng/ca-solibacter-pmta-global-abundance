#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BATCH_FILE="${1:-}"
THREADS="${THREADS:-8}"

SRA_TOOLKIT_DIR="/root/autodl-tmp/tools/sratoolkit.3.4.1-ubuntu64/bin"
if [[ -d "${SRA_TOOLKIT_DIR}" ]]; then
  export PATH="${SRA_TOOLKIT_DIR}:${PATH}"
  hash -r 2>/dev/null || true
fi

if [[ -z "${BATCH_FILE}" ]]; then
  echo "Usage: THREADS=8 bash scripts/20_process_sra_batch.sh config/sra_batches/batch_01.txt" >&2
  exit 1
fi

if [[ "${BATCH_FILE}" != /* ]]; then
  BATCH_FILE="${ROOT}/${BATCH_FILE}"
fi

if [[ ! -s "${BATCH_FILE}" ]]; then
  echo "ERROR: missing or empty batch file: ${BATCH_FILE}" >&2
  exit 1
fi

READ_DIR="${ROOT}/data/raw/reads"
MAG_DIR="${ROOT}/data/processed/solibacter_mags"
COVERM_DIR="${ROOT}/results/coverm_solibacter_all28"
PMTA_MARKER_DIR="${ROOT}/data/processed/pmta_markers"
PMTA_MARKERS="${PMTA_MARKER_DIR}/functional_pmta_marker_candidates.faa"
PMTA_DB_PREFIX="${PMTA_MARKER_DIR}/functional_pmta_marker_candidates"
PMTA_OUT_DIR="${ROOT}/results/pmta_marker_blastx_all28"
LOG_DIR="${ROOT}/logs/batches"

mkdir -p "${READ_DIR}" "${COVERM_DIR}" "${PMTA_OUT_DIR}" "${LOG_DIR}"

require_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: missing required tool: $1" >&2
    exit 1
  fi
}

require_tool prefetch
require_tool fasterq-dump
require_tool coverm
require_tool diamond
require_tool samtools
require_tool minimap2

if ! prefetch --version 2>&1 | grep -q "3.4.1"; then
  echo "WARNING: expected prefetch 3.4.1 from ${SRA_TOOLKIT_DIR}, but found:" >&2
  prefetch --version >&2 || true
  echo "If SRA download fails, run:" >&2
  echo "  export PATH=${SRA_TOOLKIT_DIR}:\$PATH" >&2
  echo "  hash -r" >&2
fi

if [[ ! -d "${MAG_DIR}" ]] || [[ -z "$(find "${MAG_DIR}" -type f -name '*.fna' -print -quit)" ]]; then
  echo "ERROR: missing Solibacter MAG FASTA files in ${MAG_DIR}" >&2
  exit 1
fi

if [[ ! -s "${PMTA_MARKERS}" ]]; then
  echo "ERROR: missing PmtA marker FASTA: ${PMTA_MARKERS}" >&2
  echo "Run first: THREADS=${THREADS} bash scripts/15_build_functional_pmta_markers.sh" >&2
  exit 1
fi

if [[ ! -s "${PMTA_DB_PREFIX}.dmnd" ]]; then
  diamond makedb --in "${PMTA_MARKERS}" -d "${PMTA_DB_PREFIX}" >/dev/null
fi

batch_name="$(basename "${BATCH_FILE}" .txt)"
echo "Processing ${batch_name} with THREADS=${THREADS}"

while read -r srr; do
  [[ -z "${srr}" ]] && continue
  echo "== ${srr} =="

  r1="${READ_DIR}/${srr}_1.fastq.gz"
  r2="${READ_DIR}/${srr}_2.fastq.gz"
  if [[ -s "${r1}" && -s "${r2}" ]]; then
    echo "Reads already exist for ${srr}; skipping download/conversion."
  else
    echo "[download] ${srr}"
    prefetch "${srr}" --output-directory "${READ_DIR}" --max-size 200G
    echo "[fasterq-dump] ${srr}"
    fasterq-dump "${READ_DIR}/${srr}" --split-files --threads "${THREADS}" -O "${READ_DIR}"
    echo "[compress] ${srr}"
    if command -v pigz >/dev/null 2>&1; then
      pigz -p "${THREADS}" -f "${READ_DIR}/${srr}"_*.fastq
    else
      gzip -f "${READ_DIR}/${srr}"_*.fastq
    fi
  fi

  coverm_out="${COVERM_DIR}/${srr}.coverm.tsv"
  if [[ -s "${coverm_out}" ]]; then
    echo "CoverM output exists for ${srr}; skipping Solibacter mapping."
  else
    echo "[coverm] ${srr}"
    coverm genome \
      --genome-fasta-directory "${MAG_DIR}" \
      --coupled "${r1}" "${r2}" \
      --methods relative_abundance trimmed_mean covered_fraction count \
      --min-covered-fraction 0 \
      --min-read-percent-identity 95 \
      --min-read-aligned-percent 90 \
      --threads "${THREADS}" \
      --output-file "${coverm_out}"
  fi

  for mate in 1 2; do
    read_file="${READ_DIR}/${srr}_${mate}.fastq.gz"
    blastx_out="${PMTA_OUT_DIR}/${srr}_${mate}.pmta_marker.blastx.tsv"
    if [[ -s "${blastx_out}" ]]; then
      echo "PmtA blastx output exists for ${srr}_${mate}; skipping."
    else
      echo "[pmtA blastx] ${srr}_${mate}"
      diamond blastx \
        -q "${read_file}" \
        -d "${PMTA_DB_PREFIX}" \
        -o "${blastx_out}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --evalue 1e-5 \
        --max-target-seqs 5 \
        --threads "${THREADS}" \
        --more-sensitive
    fi
  done
done < "${BATCH_FILE}" 2>&1 | tee "${LOG_DIR}/${batch_name}.log"

echo "Batch complete: ${batch_name}"
echo "CoverM outputs: ${COVERM_DIR}"
echo "PmtA blastx outputs: ${PMTA_OUT_DIR}"
