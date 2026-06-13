#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-8}"

MAG_DIR="${ROOT}/data/processed/solibacter_mags"
MARKER_DIR="${ROOT}/data/processed/pmta_markers"
PRODIGAL_DIR="${ROOT}/data/processed/prodigal_solibacter"
RESULT_DIR="${ROOT}/results/pmta_functional_markers"

mkdir -p "${MARKER_DIR}" "${PRODIGAL_DIR}" "${RESULT_DIR}" "${ROOT}/logs"

require_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "ERROR: missing required tool: $1" >&2
    echo "Install it in the active conda environment, e.g.:" >&2
    echo "  conda install -c bioconda -c conda-forge $1 -y" >&2
    exit 1
  fi
}

require_tool prodigal
require_tool diamond
require_tool python3

if [[ ! -d "${MAG_DIR}" ]] || [[ -z "$(find "${MAG_DIR}" -type f -name '*.fna' -print -quit)" ]]; then
  echo "ERROR: missing Solibacter MAG FASTA files in ${MAG_DIR}" >&2
  exit 1
fi

cat > "${MARKER_DIR}/Ellin6076_gene101949_functional_PmtA.faa" <<'EOF'
>Ellin6076_gene101949_functional_PmtA
MNDFMNPAEFANIRNSEEHFWWYRGMRSILFSMMERRLAGRKIDRALEAGCGTGYLSHVLQRDRGWPIVPMDLSGHGLRYARELGVQRPIQGDIRALPFAESAFDLVMSIDVIAHLPLGEEVGAARELARVTRRGGLLVVRTAALDILRSRHSEFAHERQRFTRRRLMGLFAGAGIRVLDCTYVNSLLLPAALLKFRVWEPLTGQPPESGVHPVEPWLDRLLYAPLAMEASWVGGGRSFPIGQSLVFIGERML
EOF

cat > "${MARKER_DIR}/Ellin6076_gene102575_PmtA_like_paralog.faa" <<'EOF'
>Ellin6076_gene102575_PmtA_like_paralog
MLRWIENKRLGSHLRGRGVEIGALWRRFPVAASSTVWYIDRHHSSELQREYPELKRQPVSPDVLADAGQLPFGDGSLDFVIASHILEHLPFPLACLRSWYRVLGAGGVLILKIPDMRYTFDANRRRSTLRHLVEEHLHPEEFDKRAHFQDWVEHVEGRPPGSEALRHETNRLMDLDYSIHYHVWTDEDIRELIDYTRTAMGLHWRPVLFLRAHFYRKECAVALRRE
EOF

echo "[1/4] Predicting proteins from Solibacter MAGs with Prodigal..."
for f in "${MAG_DIR}"/*.fna; do
  base="$(basename "${f}" .fna)"
  raw_faa="${PRODIGAL_DIR}/${base}.raw.faa"
  prefixed_faa="${PRODIGAL_DIR}/${base}.faa"
  ffn="${PRODIGAL_DIR}/${base}.ffn"
  if [[ ! -s "${prefixed_faa}" ]]; then
    prodigal -i "${f}" -a "${raw_faa}" -d "${ffn}" -p meta -q
    awk -v b="${base}" '/^>/{sub(/^>/, ">" b "|")} {print}' "${raw_faa}" > "${prefixed_faa}"
  fi
done

cat "${PRODIGAL_DIR}"/*.faa > "${MARKER_DIR}/solibacter_predicted_proteins.faa"

echo "[2/4] Building DIAMOND database for Solibacter predicted proteins..."
diamond makedb \
  --in "${MARKER_DIR}/solibacter_predicted_proteins.faa" \
  -d "${MARKER_DIR}/solibacter_predicted_proteins.dmnd" \
  >/dev/null

echo "[3/4] Searching Solibacter proteins with functional gene101949 PmtA seed..."
diamond blastp \
  -q "${MARKER_DIR}/Ellin6076_gene101949_functional_PmtA.faa" \
  -d "${MARKER_DIR}/solibacter_predicted_proteins.dmnd" \
  -o "${RESULT_DIR}/gene101949_vs_solibacter_proteins.tsv" \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
  --evalue 1e-10 \
  --max-target-seqs 10000 \
  --threads "${THREADS}" \
  --sensitive

echo "[4/4] Filtering hits by protein-level PmtA features..."
python3 "${ROOT}/scripts/16_filter_pmta_protein_hits.py" \
  --blast "${RESULT_DIR}/gene101949_vs_solibacter_proteins.tsv" \
  --proteins "${MARKER_DIR}/solibacter_predicted_proteins.faa" \
  --functional-seed "${MARKER_DIR}/Ellin6076_gene101949_functional_PmtA.faa" \
  --paralog-seed "${MARKER_DIR}/Ellin6076_gene102575_PmtA_like_paralog.faa" \
  --out-tsv "${RESULT_DIR}/functional_pmta_marker_candidates.tsv" \
  --out-faa "${MARKER_DIR}/functional_pmta_marker_candidates.faa"

echo "Done."
echo "Candidate table: ${RESULT_DIR}/functional_pmta_marker_candidates.tsv"
echo "Marker FASTA:    ${MARKER_DIR}/functional_pmta_marker_candidates.faa"
