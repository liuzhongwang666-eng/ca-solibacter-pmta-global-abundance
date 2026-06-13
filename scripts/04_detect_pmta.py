from pathlib import Path
import argparse

ROOT = Path(__file__).resolve().parents[1]
PROTEINS = ROOT / "data" / "processed" / "solibacter_proteins.faa"
PMTA = ROOT / "data" / "processed" / "Ellin6076_pmtA_seed.faa"
OUT = ROOT / "results" / "pmta_detection_commands.sh"

COMMANDS = f"""#!/usr/bin/env bash
set -euo pipefail

# 1. Predict proteins from selected Solibacter MAGs if proteins are not available.
# Example:
# mkdir -p data/processed/prodigal
# for f in data/processed/solibacter_mags/*.fa*; do
#   b=$(basename "$f")
#   prodigal -i "$f" -a "data/processed/prodigal/${{b}}.faa" -d "data/processed/prodigal/${{b}}.ffn" -p meta
# done
# cat data/processed/prodigal/*.faa > data/processed/solibacter_proteins.faa

# 2. DIAMOND/BLASTP search using Ellin6076 PmtA as seed.
diamond makedb --in {PROTEINS.as_posix()} -d data/processed/solibacter_proteins.dmnd
diamond blastp \\
  -q {PMTA.as_posix()} \\
  -d data/processed/solibacter_proteins.dmnd \\
  -o results/pmta_diamond_hits.tsv \\
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \\
  --evalue 1e-20 --max-target-seqs 10000

# 3. Optional HMMER confirmation after building a PmtA HMM from curated homologs.
# hmmbuild data/processed/pmtA_curated.hmm data/processed/pmtA_curated_alignment.sto
# hmmsearch --tblout results/pmta_hmmer_hits.tbl --cpu 8 data/processed/pmtA_curated.hmm {PROTEINS.as_posix()}

# 4. Filter hits downstream:
# identity >= 40%, query coverage >= 70%, evalue <= 1e-20.
"""


parser = argparse.ArgumentParser()
parser.add_argument("--write-commands", action="store_true")
args = parser.parse_args()

OUT.write_text(COMMANDS, encoding="utf-8")
print(f"Wrote pmtA detection command template to {OUT}")
if not args.write_commands:
    print("This script writes commands only. Run the generated shell script in a Unix/HPC environment.")
