#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${1:-solibacter_mapping}"

if command -v mamba >/dev/null 2>&1; then
  mamba create -n "${ENV_NAME}" -c bioconda -c conda-forge \
    coverm bwa minimap2 samtools sra-tools seqkit pigz curl wget diamond prodigal python=3.11 -y
elif command -v conda >/dev/null 2>&1; then
  conda create -n "${ENV_NAME}" -c bioconda -c conda-forge \
    coverm bwa minimap2 samtools sra-tools seqkit pigz curl wget diamond prodigal python=3.11 -y
else
  echo "ERROR: conda or mamba is required on the HPC/server." >&2
  exit 1
fi

echo "Activate with: conda activate ${ENV_NAME}"
