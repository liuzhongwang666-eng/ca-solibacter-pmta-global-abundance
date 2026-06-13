#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="${ROOT}/results/candidate_solibacter_manifest.tsv"
MAG_ROOT="${ROOT}/data/raw/mag"
OUT_DIR="${ROOT}/data/processed/solibacter_mags"
ID_LIST="${ROOT}/data/processed/solibacter_mag_ids.txt"
LOG="${ROOT}/logs/extract_solibacter_mags.log"

mkdir -p "${OUT_DIR}" "$(dirname "${ID_LIST}")" "$(dirname "${LOG}")"

if [[ ! -s "${MANIFEST}" ]]; then
  echo "ERROR: Missing ${MANIFEST}" >&2
  exit 1
fi
if [[ ! -d "${MAG_ROOT}" ]]; then
  echo "ERROR: Missing MAG directory ${MAG_ROOT}. Run 08_download_smag_mag_archive.sh first." >&2
  exit 1
fi

awk -F'\t' 'NR>1 {print $2}' "${MANIFEST}" | sort -u > "${ID_LIST}"
: > "${LOG}"

while read -r mag_id; do
  [[ -z "${mag_id}" ]] && continue
  match="$(find "${MAG_ROOT}" -type f \( -name "*${mag_id}*.fa" -o -name "*${mag_id}*.fna" -o -name "*${mag_id}*.fasta" -o -name "*${mag_id}*.fa.gz" -o -name "*${mag_id}*.fna.gz" -o -name "*${mag_id}*.fasta.gz" \) | head -n 1 || true)"
  if [[ -z "${match}" ]]; then
    echo -e "${mag_id}\tMISSING" | tee -a "${LOG}"
    continue
  fi
  if [[ "${match}" == *.gz ]]; then
    gzip -cd "${match}" > "${OUT_DIR}/${mag_id}.fna"
  else
    cp -f "${match}" "${OUT_DIR}/${mag_id}.fna"
  fi
  echo -e "${mag_id}\tCOPIED\t${match}" | tee -a "${LOG}"
done < "${ID_LIST}"

echo "Extracted $(find "${OUT_DIR}" -type f -name '*.fna' | wc -l) Solibacter MAG FASTA files to ${OUT_DIR}"
