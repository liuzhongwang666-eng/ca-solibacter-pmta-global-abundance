from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "results" / "solibacter_sample_manifest.tsv"
OUT = ROOT / "results" / "sra_accessions_for_solibacter_candidates.txt"

with META.open(encoding="utf-8", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

accessions = sorted({row["sample"] for row in rows if row.get("sample")})
OUT.write_text("\n".join(accessions) + "\n", encoding="utf-8")
print(f"Wrote {len(accessions)} SRA accessions to {OUT}")
