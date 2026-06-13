from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
CANDIDATES = ROOT / "results" / "candidate_solibacter_mags.tsv"
MAG_DIR = ROOT / "data" / "raw" / "mag"
OUT_DIR = ROOT / "data" / "processed" / "solibacter_mags"
MANIFEST = ROOT / "data" / "processed" / "solibacter_mag_manifest.tsv"


def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


try:
    import pandas as pd
except Exception:
    fail("This script requires pandas.")

if not CANDIDATES.exists():
    fail("Run scripts/02_screen_smag_taxonomy.py first.")
if not MAG_DIR.exists():
    fail("Extract SMAG MAG FASTA files into data/raw/mag/ first.")

df = pd.read_csv(CANDIDATES, sep="\t")
if df.empty:
    fail("Candidate table is empty.")

candidate_text = df.astype(str)
ids = set()
for col in candidate_text.columns:
    lower = col.lower()
    if any(k in lower for k in ["mag", "genome", "bin", "sgb", "id"]):
        ids.update(v for v in candidate_text[col].dropna().astype(str) if v and v.lower() != "nan")

if not ids:
    fail("Could not infer MAG IDs. Rename the MAG ID column or manually create a manifest.")

OUT_DIR.mkdir(parents=True, exist_ok=True)
rows = ["candidate_id\tmatched_file\tstatus"]
for cid in sorted(ids):
    matches = list(MAG_DIR.rglob(f"*{cid}*"))
    fasta_matches = [p for p in matches if p.suffix.lower() in {".fa", ".fna", ".fasta"} or ".fa" in p.name]
    if not fasta_matches:
        rows.append(f"{cid}\tNA\tmissing")
        continue
    src = fasta_matches[0]
    dest = OUT_DIR / src.name
    if not dest.exists():
        shutil.copy2(src, dest)
    rows.append(f"{cid}\t{src}\tcopied")

MANIFEST.write_text("\n".join(rows) + "\n", encoding="utf-8")
print(f"Wrote manifest to {MANIFEST}")
