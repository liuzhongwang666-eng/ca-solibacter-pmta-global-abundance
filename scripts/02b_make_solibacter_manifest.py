from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "candidate_solibacter_mags.tsv"
OUT = ROOT / "results" / "candidate_solibacter_manifest.tsv"

FIELDS = {
    "smag_index": "Number",
    "mag_id": "user_genome",
    "completeness": "completeness",
    "contamination": "contamination",
    "genome_length_bp": "length",
    "n50": "N50",
    "gtdb_domain": "GTDB202_classification",
    "gtdb_phylum": "GTDB202_phylum",
    "gtdb_class": "GTDB202_class",
    "gtdb_order": "GTDB202_order",
    "gtdb_family": "GTDB202_family",
    "gtdb_genus": "GTDB202_genus",
    "gtdb_species": "GTDB202_species",
    "closest_references": "other_related_references.genome_id.species_name.radius.ANI.AF.",
    "derep_sgb_taxonomy": "GTDB214",
    "checkm2_completeness": "checkm2_completeness",
    "checkm2_contamination": "checkm2_contamination",
    "genome_size_bp_checkm2": "Genome_Size",
    "total_coding_sequences": "Total_Coding_Sequences",
}


def qc_label(completeness, contamination):
    try:
        comp = float(completeness)
        cont = float(contamination)
    except Exception:
        return "unknown"
    if comp >= 90 and cont < 5:
        return "high_quality"
    if comp >= 70 and cont < 10:
        return "medium_quality"
    return "low_or_partial_quality"


with IN.open(encoding="utf-8", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

out_fields = list(FIELDS.keys()) + ["quality_group", "source_file", "source_sheet", "row_number"]
with OUT.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=out_fields, delimiter="\t")
    writer.writeheader()
    for row in rows:
        out = {name: row.get(col, "") for name, col in FIELDS.items()}
        out["quality_group"] = qc_label(out["completeness"], out["contamination"])
        out["source_file"] = row.get("source_file", "")
        out["source_sheet"] = row.get("source_sheet", "")
        out["row_number"] = row.get("row_number", "")
        writer.writerow(out)

print(f"Wrote {len(rows)} rows to {OUT}")
