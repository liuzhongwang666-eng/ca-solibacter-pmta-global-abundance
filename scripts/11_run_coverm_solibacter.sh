#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:-pilot}"

if [[ "${MODE}" == "pilot" ]]; then
  ACCESSIONS="${ROOT}/config/pilot_sra_accessions.txt"
elif [[ "${MODE}" == "all28" ]]; then
  ACCESSIONS="${ROOT}/results/sra_accessions_for_solibacter_candidates.txt"
else
  echo "Usage: bash scripts/11_run_coverm_solibacter.sh [pilot|all28]" >&2
  exit 1
fi

MAG_DIR="${ROOT}/data/processed/solibacter_mags"
READ_DIR="${ROOT}/data/raw/reads"
OUT_DIR="${ROOT}/results/coverm_solibacter_${MODE}"
mkdir -p "${OUT_DIR}" "${ROOT}/logs"

if [[ ! -d "${MAG_DIR}" ]] || [[ -z "$(find "${MAG_DIR}" -type f -name '*.fna' -print -quit)" ]]; then
  echo "ERROR: Missing Solibacter MAG FASTA files in ${MAG_DIR}" >&2
  exit 1
fi

while read -r srr; do
  [[ -z "${srr}" ]] && continue
  r1="${READ_DIR}/${srr}_1.fastq.gz"
  r2="${READ_DIR}/${srr}_2.fastq.gz"
  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "ERROR: Missing paired reads for ${srr}: ${r1} ${r2}" >&2
    exit 1
  fi
  echo "Running CoverM for ${srr}"
  coverm genome \
    --genome-fasta-directory "${MAG_DIR}" \
    --coupled "${r1}" "${r2}" \
    --methods relative_abundance trimmed_mean covered_fraction count \
    --min-covered-fraction 0 \
    --min-read-percent-identity 95 \
    --min-read-aligned-percent 90 \
    --threads "${THREADS:-8}" \
    --output-file "${OUT_DIR}/${srr}.coverm.tsv"
done < "${ACCESSIONS}"

python3 "${ROOT}/scripts/12_summarize_solibacter_abundance.py" \
  --coverm-dir "${OUT_DIR}" \
  --metadata "${ROOT}/results/solibacter_sample_manifest.tsv" \
  --out "${ROOT}/results/sample_solibacter_abundance_${MODE}.tsv"

echo "Sample-level abundance written to results/sample_solibacter_abundance_${MODE}.tsv"
