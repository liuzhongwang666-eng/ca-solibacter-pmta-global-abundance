#!/usr/bin/env python3
"""Fetch SRA runs for soil-confirmed NCBI Candidatus Solibacter source projects.

NCBI MAGs are used only as project-discovery clues. This script retrieves raw
metagenomic runs from their source BioProjects and flags whether a run matches a
BioSample associated with a Solibacter assembly.
"""

import argparse
import csv
import io
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path


RUNINFO_URL = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo"


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fetch_runinfo(project, timeout=120):
    query = urllib.parse.urlencode({"acc": project})
    with urllib.request.urlopen(f"{RUNINFO_URL}?{query}", timeout=timeout) as resp:
        text = resp.read().decode("utf-8", errors="replace")
    if not text.strip() or text.startswith("Error"):
        return []
    return list(csv.DictReader(io.StringIO(text)))


def is_metagenomic(row):
    strategy = (row.get("LibraryStrategy") or "").upper()
    source = (row.get("LibrarySource") or "").upper()
    return strategy in {"WGS", "METAGENOMIC", "AMPLICON"} or "METAGENOMIC" in source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--projects", default="config/ncbi_solibacter_soil_source_projects.tsv")
    parser.add_argument("--out", default="results/ncbi_solibacter_soil_project_sra_runs.tsv")
    parser.add_argument("--excluded-out", default="results/ncbi_solibacter_soil_project_sra_runs_excluded.tsv")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between NCBI requests.")
    args = parser.parse_args()

    projects = read_tsv(args.projects)
    source_biosamples = {}
    project_meta = {}
    for row in projects:
        project = row["bioproject"].strip()
        source_biosamples.setdefault(project, set()).add(row.get("biosample", "").strip())
        project_meta.setdefault(project, []).append(row)

    included = []
    excluded = []
    for project in sorted(project_meta):
        print(f"[NCBI runinfo] {project}", file=sys.stderr)
        try:
            runs = fetch_runinfo(project)
        except Exception as exc:  # noqa: BLE001
            excluded.append(
                {
                    "bioproject": project,
                    "run": "",
                    "biosample": "",
                    "reason": f"runinfo_fetch_failed:{exc}",
                }
            )
            continue
        time.sleep(args.sleep)
        wanted_biosamples = source_biosamples.get(project, set())
        for run in runs:
            run_acc = run.get("Run", "")
            biosample = run.get("BioSample", "")
            strategy = run.get("LibraryStrategy", "")
            source = run.get("LibrarySource", "")
            if not run_acc:
                continue
            matched_source_biosample = biosample in wanted_biosamples
            if not is_metagenomic(run):
                excluded.append(
                    {
                        "bioproject": project,
                        "run": run_acc,
                        "biosample": biosample,
                        "reason": f"non_metagenomic_strategy:{strategy}/{source}",
                    }
                )
                continue
            template = project_meta[project][0]
            included.append(
                {
                    "run": run_acc,
                    "bioproject": project,
                    "biosample": biosample,
                    "matched_ncbi_solibacter_biosample": "yes" if matched_source_biosample else "no",
                    "source_type": template.get("source_type", "NCBI_Solibacter_source_project_reads"),
                    "ecosystem": template.get("ecosystem", ""),
                    "location": template.get("location", ""),
                    "latitude": template.get("latitude", ""),
                    "longitude": template.get("longitude", ""),
                    "soil_only": "yes",
                    "acidic_soil": template.get("acidic_soil", ""),
                    "soil_source_evidence": template.get("soil_source_evidence", ""),
                    "pH": template.get("pH", ""),
                    "pH_source": template.get("pH_source", ""),
                    "recommendation": template.get("recommendation", ""),
                    "overlap_class": template.get("overlap_class", ""),
                    "LibraryStrategy": strategy,
                    "LibrarySource": source,
                    "spots": run.get("spots", ""),
                    "bases": run.get("bases", ""),
                    "size_MB": run.get("size_MB", ""),
                    "AssemblyName": run.get("AssemblyName", ""),
                }
            )

    for path, rows in ((args.out, included), (args.excluded_out, excluded)):
        out = Path(path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fields = sorted({k for row in rows for k in row}) if rows else ["run", "bioproject", "biosample", "reason"]
        with out.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {out}")


if __name__ == "__main__":
    main()
