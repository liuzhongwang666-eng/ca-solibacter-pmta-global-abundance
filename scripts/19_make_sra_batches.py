#!/usr/bin/env python3
import argparse
from pathlib import Path


parser = argparse.ArgumentParser(description="Split SRA accessions into fixed-size batch files.")
parser.add_argument(
    "--accessions",
    default="results/sra_accessions_for_solibacter_candidates.txt",
    help="Input accession list, one SRR per line.",
)
parser.add_argument(
    "--out-dir",
    default="config/sra_batches",
    help="Output directory for batch accession files.",
)
parser.add_argument("--batch-size", type=int, default=5)
parser.add_argument(
    "--skip",
    default="config/pilot_sra_accessions.txt",
    help="Optional accession list to skip, e.g. already completed pilot samples. Use empty string to disable.",
)
args = parser.parse_args()

root = Path.cwd()
accession_path = root / args.accessions
out_dir = root / args.out_dir

accessions = [
    line.strip()
    for line in accession_path.read_text(encoding="utf-8").splitlines()
    if line.strip()
]

skip = set()
if args.skip:
    skip_path = root / args.skip
    if skip_path.exists():
        skip = {
            line.strip()
            for line in skip_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        }

remaining = [srr for srr in accessions if srr not in skip]
out_dir.mkdir(parents=True, exist_ok=True)

for old in out_dir.glob("batch_*.txt"):
    old.unlink()

for i in range(0, len(remaining), args.batch_size):
    batch = remaining[i : i + args.batch_size]
    batch_no = i // args.batch_size + 1
    path = out_dir / f"batch_{batch_no:02d}.txt"
    path.write_text("\n".join(batch) + "\n", encoding="utf-8")
    print(f"Wrote {path} ({len(batch)} samples)")

print(f"Input accessions: {len(accessions)}")
print(f"Skipped accessions: {len(skip)}")
print(f"Remaining accessions: {len(remaining)}")
