#!/usr/bin/env bash
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT}" || exit 1

echo "== Remote pilot preflight check =="
echo "Project root: ${ROOT}"
echo

echo "== Host and working directory =="
hostname || true
pwd
echo

echo "== Disk space =="
df -h .
echo

echo "== Pilot accessions =="
if [[ -s config/pilot_sra_accessions.txt ]]; then
  cat config/pilot_sra_accessions.txt
else
  echo "MISSING: config/pilot_sra_accessions.txt"
fi
echo

echo "== Required metadata/results files =="
for f in \
  results/candidate_solibacter_manifest.tsv \
  results/solibacter_sample_manifest.tsv \
  results/sra_accessions_for_solibacter_candidates.txt
do
  if [[ -s "$f" ]]; then
    echo "OK: $f"
  else
    echo "MISSING: $f"
  fi
done
echo

echo "== Conda/mamba =="
if command -v mamba >/dev/null 2>&1; then
  echo "OK: mamba $(command -v mamba)"
elif command -v conda >/dev/null 2>&1; then
  echo "OK: conda $(command -v conda)"
else
  echo "WARNING: conda/mamba not found. Try:"
  echo "  source ~/miniconda3/etc/profile.d/conda.sh"
  echo "  source /root/miniconda3/etc/profile.d/conda.sh"
fi
echo

echo "== Mapping tools =="
for tool in coverm bwa minimap2 samtools prefetch fasterq-dump seqkit pigz gzip; do
  if command -v "$tool" >/dev/null 2>&1; then
    echo "OK: $tool"
  else
    echo "NOT_FOUND: $tool"
  fi
done
echo

echo "== Existing MAG/reads status =="
if [[ -d data/raw/mag ]]; then
  echo "data/raw/mag exists: $(find data/raw/mag -type f | wc -l) files"
else
  echo "data/raw/mag not found yet. Run scripts/08_download_smag_mag_archive.sh"
fi

if [[ -d data/processed/solibacter_mags ]]; then
  echo "data/processed/solibacter_mags exists: $(find data/processed/solibacter_mags -type f -name '*.fna' | wc -l) .fna files"
else
  echo "data/processed/solibacter_mags not found yet. Run scripts/09_extract_solibacter_mags_hpc.sh after MAG download."
fi

if [[ -d data/raw/reads ]]; then
  echo "data/raw/reads exists: $(find data/raw/reads -type f -name '*.fastq.gz' | wc -l) fastq.gz files"
else
  echo "data/raw/reads not found yet. Run scripts/10_download_sra_reads.sh pilot"
fi
echo

echo "== Recommended next command =="
echo "If conda/mamba and tools are missing:"
echo "  bash scripts/07_install_mapping_env.sh"
echo "  conda activate solibacter_mapping"
echo
echo "Then run:"
echo "  bash scripts/08_download_smag_mag_archive.sh"
echo "  bash scripts/09_extract_solibacter_mags_hpc.sh"
echo "  THREADS=8 bash scripts/10_download_sra_reads.sh pilot"
echo "  THREADS=8 bash scripts/11_run_coverm_solibacter.sh pilot"
