from pathlib import Path
import argparse
import csv
import re


def clean_name(name):
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_").lower()


def read_table(path):
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def find_column(fieldnames, keyword):
    for field in fieldnames:
        if keyword in clean_name(field):
            return field
    return None


parser = argparse.ArgumentParser()
parser.add_argument("--coverm-dir", required=True)
parser.add_argument("--metadata", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

coverm_dir = Path(args.coverm_dir)
metadata_path = Path(args.metadata)
out_path = Path(args.out)

metadata_rows = read_table(metadata_path)
metadata_by_sample = {row.get("sample", ""): row for row in metadata_rows}

rows_out = []
for path in sorted(coverm_dir.glob("*.coverm.tsv")):
    sample = path.name.split(".coverm.tsv")[0]
    rows = read_table(path)
    if not rows:
        continue
    fields = rows[0].keys()
    rel_col = find_column(fields, "relative_abundance")
    trim_col = find_column(fields, "trimmed_mean")
    covered_col = find_column(fields, "covered_fraction")
    count_col = find_column(fields, "count")

    def sum_float(col):
        if not col:
            return ""
        total = 0.0
        for row in rows:
            try:
                total += float(row.get(col, "0") or 0)
            except ValueError:
                pass
        return total

    def mean_float(col):
        if not col:
            return ""
        vals = []
        for row in rows:
            try:
                vals.append(float(row.get(col, "0") or 0))
            except ValueError:
                pass
        return sum(vals) / len(vals) if vals else ""

    meta = metadata_by_sample.get(sample, {})
    rows_out.append(
        {
            "sample_id": sample,
            "Candidatus_Solibacter_relative_abundance_sum": sum_float(rel_col),
            "Candidatus_Solibacter_trimmed_mean_sum": sum_float(trim_col),
            "Candidatus_Solibacter_covered_fraction_mean": mean_float(covered_col),
            "Candidatus_Solibacter_mapped_count_sum": sum_float(count_col),
            "n_solibacter_genomes": len(rows),
            "latitude": meta.get("latitude", ""),
            "longitude": meta.get("longitude", ""),
            "ecosystem": meta.get("ecosystem", ""),
            "continent": meta.get("continent", ""),
            "nation": meta.get("nation", ""),
            "bioproject": meta.get("bioproject", ""),
            "biosample": meta.get("biosample", ""),
        }
    )

out_path.parent.mkdir(parents=True, exist_ok=True)
fieldnames = [
    "sample_id",
    "Candidatus_Solibacter_relative_abundance_sum",
    "Candidatus_Solibacter_trimmed_mean_sum",
    "Candidatus_Solibacter_covered_fraction_mean",
    "Candidatus_Solibacter_mapped_count_sum",
    "n_solibacter_genomes",
    "latitude",
    "longitude",
    "ecosystem",
    "continent",
    "nation",
    "bioproject",
    "biosample",
]
with out_path.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows_out)

print(f"Wrote {len(rows_out)} sample abundance rows to {out_path}")
