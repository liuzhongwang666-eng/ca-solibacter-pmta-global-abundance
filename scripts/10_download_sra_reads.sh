#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:-pilot}"

if [[ "${MODE}" == "pilot" ]]; then
  ACCESSIONS="${ROOT}/config/pilot_sra_accessions.txt"
elif [[ "${MODE}" == "all28" ]]; then
  ACCESSIONS="${ROOT}/results/sra_accessions_for_solibacter_candidates.txt"
else
  echo "Usage: bash scripts/10_download_sra_reads.sh [pilot|all28]" >&2
  exit 1
fi

READ_DIR="${ROOT}/data/raw/reads"
mkdir -p "${READ_DIR}" "${ROOT}/logs"

while read -r srr; do
  [[ -z "${srr}" ]] && continue
  echo "Downloading ${srr}"
  if prefetch "${srr}" --output-directory "${READ_DIR}"; then
    :
  else
    echo "prefetch --output-directory failed; retrying with legacy -O syntax for ${srr}"
    prefetch "${srr}" -O "${READ_DIR}"
  fi
  fasterq-dump "${READ_DIR}/${srr}" --split-files --threads "${THREADS:-8}" -O "${READ_DIR}"
  if command -v pigz >/dev/null 2>&1; then
    pigz -p "${THREADS:-8}" -f "${READ_DIR}/${srr}"_*.fastq
  else
    gzip -f "${READ_DIR}/${srr}"_*.fastq
  fi
done < "${ACCESSIONS}"

echo "Reads downloaded to ${READ_DIR}"
