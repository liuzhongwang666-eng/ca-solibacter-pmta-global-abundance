#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MODE="${1:-pilot}"
THREADS="${THREADS:-8}"

if [[ "${MODE}" == "pilot" ]]; then
  ACCESSIONS="${ROOT}/config/pilot_sra_accessions.txt"
elif [[ "${MODE}" == "all28" ]]; then
  ACCESSIONS="${ROOT}/results/sra_accessions_for_solibacter_candidates.txt"
else
  echo "Usage: bash scripts/17_run_pmta_marker_blastx_pilot.sh [pilot|all28]" >&2
  exit 1
fi

READ_DIR="${ROOT}/data/raw/reads"
MARKER_DIR="${ROOT}/data/processed/pmta_markers"
MARKERS="${MARKER_DIR}/functional_pmta_marker_candidates.faa"
DB_PREFIX="${MARKER_DIR}/functional_pmta_marker_candidates"
OUT_DIR="${ROOT}/results/pmta_marker_blastx_${MODE}"

mkdir -p "${OUT_DIR}" "${ROOT}/logs"

if ! command -v diamond >/dev/null 2>&1; then
  echo "ERROR: missing required tool: diamond" >&2
  echo "Install it in the active conda environment:" >&2
  echo "  conda install -c bioconda -c conda-forge diamond -y" >&2
  exit 1
fi

if [[ ! -s "${MARKERS}" ]]; then
  echo "ERROR: missing PmtA marker FASTA: ${MARKERS}" >&2
  echo "Run first: THREADS=${THREADS} bash scripts/15_build_functional_pmta_markers.sh" >&2
  exit 1
fi

echo "[1/3] Building DIAMOND database for functional PmtA markers..."
diamond makedb --in "${MARKERS}" -d "${DB_PREFIX}" >/dev/null

echo "[2/3] Running translated read searches against functional PmtA markers..."
while read -r srr; do
  [[ -z "${srr}" ]] && continue
  r1="${READ_DIR}/${srr}_1.fastq.gz"
  r2="${READ_DIR}/${srr}_2.fastq.gz"
  if [[ ! -s "${r1}" || ! -s "${r2}" ]]; then
    echo "ERROR: missing paired reads for ${srr}: ${r1} ${r2}" >&2
    exit 1
  fi
  for mate in 1 2; do
    read_file="${READ_DIR}/${srr}_${mate}.fastq.gz"
    out_file="${OUT_DIR}/${srr}_${mate}.pmta_marker.blastx.tsv"
    echo "Running DIAMOND blastx for ${srr}_${mate}"
    diamond blastx \
      -q "${read_file}" \
      -d "${DB_PREFIX}" \
      -o "${out_file}" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --evalue 1e-5 \
      --max-target-seqs 5 \
      --threads "${THREADS}" \
      --more-sensitive
  done
done < "${ACCESSIONS}"

echo "[3/3] Summarizing PmtA marker abundance..."
python3 "${ROOT}/scripts/18_summarize_pmta_marker_abundance.py" \
  --blastx-dir "${OUT_DIR}" \
  --reads-dir "${READ_DIR}" \
  --samples "${ACCESSIONS}" \
  --markers "${MARKERS}" \
  --out "${ROOT}/results/sample_pmta_marker_abundance_${MODE}.csv"

echo "Sample-level PmtA marker abundance written to results/sample_pmta_marker_abundance_${MODE}.csv"
