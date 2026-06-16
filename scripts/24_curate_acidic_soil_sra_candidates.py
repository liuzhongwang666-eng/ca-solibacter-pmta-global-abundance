#!/usr/bin/env python3
"""Curate soil-only acidic metagenome SRA candidates for extension mapping."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_table(path, delimiter="\t"):
    p = Path(path)
    if not p.exists():
        return []
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def truthy(value):
    return str(value or "").strip().lower() in {"yes", "true", "1", "y"}


def float_or_none(value):
    try:
        if value in ("", None):
            return None
        return float(value)
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ncbi-runs", default="results/ncbi_solibacter_soil_project_sra_runs.tsv")
    parser.add_argument("--manual-projects", default="config/acidic_soil_project_candidates.tsv")
    parser.add_argument("--existing-samples", default="results/sra_accessions_for_solibacter_candidates.txt")
    parser.add_argument("--out", default="results/acidic_soil_sra_extension_candidates.tsv")
    parser.add_argument("--excluded-out", default="results/acidic_soil_sra_extension_excluded.tsv")
    parser.add_argument("--pilot-out", default="config/acidic_soil_pilot_sra_accessions.txt")
    parser.add_argument("--max-per-project", type=int, default=5)
    parser.add_argument("--max-projects", type=int, default=5)
    parser.add_argument("--ph-threshold", type=float, default=5.5)
    args = parser.parse_args()

    existing = set()
    existing_path = Path(args.existing_samples)
    if existing_path.exists():
        existing = {x.strip() for x in existing_path.read_text(encoding="utf-8").splitlines() if x.strip()}

    rows = []
    for row in read_table(args.ncbi_runs):
        row = dict(row)
        row.setdefault("source_type", "NCBI_Solibacter_source_project_reads")
        row.setdefault("curation_status", "candidate")
        row.setdefault("soil_only", "yes")
        rows.append(row)
    for row in read_table(args.manual_projects):
        row = dict(row)
        row.setdefault("source_type", "external_acidic_soil_project_reads")
        rows.append(row)

    included = []
    excluded = []
    seen = set()
    for row in rows:
        run = row.get("run", "").strip()
        if not run:
            row["exclusion_reason"] = "missing_run_accession"
            excluded.append(row)
            continue
        if run in seen:
            row["exclusion_reason"] = "duplicate_run"
            excluded.append(row)
            continue
        seen.add(run)
        if run in existing:
            row["exclusion_reason"] = "already_in_existing_28_sra"
            excluded.append(row)
            continue
        if not truthy(row.get("soil_only", "yes")):
            row["exclusion_reason"] = "not_soil_only"
            excluded.append(row)
            continue
        ph = float_or_none(row.get("pH") or row.get("soil_pH"))
        acidic_evidence = truthy(row.get("acidic_soil", ""))
        if ph is not None and ph <= args.ph_threshold:
            acidic_evidence = True
        if not acidic_evidence:
            row["exclusion_reason"] = "no_acidic_soil_evidence"
            excluded.append(row)
            continue
        row["record_type"] = row.get("source_type", "external_acidic_soil_project_reads")
        row["exclusion_reason"] = ""
        included.append(row)

    by_project = defaultdict(list)
    for row in included:
        by_project[row.get("bioproject", "unknown")].append(row)

    selected = []
    for project in sorted(by_project):
        if len(selected) >= args.max_projects * args.max_per_project:
            break
        project_rows = sorted(
            by_project[project],
            key=lambda r: (
                r.get("matched_ncbi_solibacter_biosample") != "yes",
                -(float_or_none(r.get("bases")) or 0),
                r.get("run", ""),
            ),
        )
        selected.extend(project_rows[: args.max_per_project])
        if len({r.get("bioproject") for r in selected}) >= args.max_projects:
            continue

    for path, table in ((args.out, included), (args.excluded_out, excluded)):
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fields = sorted({k for row in table for k in row}) if table else ["run", "bioproject", "exclusion_reason"]
        with out.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(table)
        print(f"Wrote {len(table)} rows to {out}")

    pilot = Path(args.pilot_out)
    pilot.parent.mkdir(parents=True, exist_ok=True)
    pilot.write_text("\n".join(row["run"] for row in selected) + ("\n" if selected else ""), encoding="ascii")
    print(f"Wrote {len(selected)} pilot accessions to {pilot}")


if __name__ == "__main__":
    main()
