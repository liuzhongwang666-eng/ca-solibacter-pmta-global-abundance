#!/usr/bin/env python3
"""Build a BioProject-deduplicated global acidic-soil pilot SRA list.

The script searches NCBI SRA by continent/acidic-soil type, retrieves
BioSample metadata, removes non-soil records, limits each BioProject to at most
two selected runs per group, and writes a mapping-ready accession list.
"""

from __future__ import annotations

import argparse
import csv
import io
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path


ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


QUERY_GROUPS = [
    {
        "group": "Asia red/acid soil",
        "continent": "Asia",
        "acidic_soil_type": "China red/acid soil",
        "ecosystem": "Red/acid soil",
        "target": 5,
        "queries": [
            '("acidic soil" OR "red soil") AND China AND (metagenomic OR metagenome)',
            '(red soil OR Acrisol OR Ferralsol) AND China AND soil AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "Asia acidic paddy soil",
        "continent": "Asia",
        "acidic_soil_type": "Acidic paddy soil",
        "ecosystem": "Paddy soil",
        "target": 5,
        "queries": [
            '("paddy soil" AND (acidic OR acid OR red soil)) AND (metagenomic OR metagenome)',
            '("rice paddy" AND acidic AND soil) AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "Europe acidic forest soil",
        "continent": "Europe",
        "acidic_soil_type": "Acidic forest soil",
        "ecosystem": "Forest soil",
        "target": 5,
        "queries": [
            '("acidic forest soil" OR "forest soil") AND Europe AND (metagenomic OR metagenome)',
            '(podzol OR acidic) AND forest soil AND Europe AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "North America peat/permafrost soil",
        "continent": "North America",
        "acidic_soil_type": "Peat/permafrost/tundra soil",
        "ecosystem": "Peat/permafrost soil",
        "target": 5,
        "queries": [
            '((peatland OR permafrost OR tundra OR bog) AND soil) AND (USA OR "United States" OR Canada OR Alaska OR Minnesota OR "North America") AND (metagenomic OR metagenome)',
            '(bog OR peat OR permafrost) AND soil AND (USA OR "United States" OR Canada) AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "South America Amazon/tropical soil",
        "continent": "South America",
        "acidic_soil_type": "Amazon/tropical acidic soil candidate",
        "ecosystem": "Amazon/tropical soil",
        "target": 5,
        "queries": [
            '(Amazon OR Amazonian OR "tropical forest") AND soil AND (metagenomic OR metagenome)',
            '(acrisol OR ferralsol OR oxisol) AND (Brazil OR Peru OR Colombia OR Ecuador) AND (metagenomic OR metagenome)',
            'Brazil AND (forest soil OR pasture soil) AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "Africa acidic/Ferralsol/Acrisol soil",
        "continent": "Africa",
        "acidic_soil_type": "Acidic/Ferralsol/Acrisol soil candidate",
        "ecosystem": "Acidic/Ferralsol/Acrisol soil candidate",
        "target": 5,
        "queries": [
            'Africa AND ("acidic soil" OR "acid soil" OR ferralsol OR acrisol) AND (metagenomic OR metagenome)',
            '(ferralsol OR acrisol OR oxisol) AND Africa AND soil AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "Oceania acid sulfate soil",
        "continent": "Oceania",
        "acidic_soil_type": "Acid sulfate soil",
        "ecosystem": "Acid sulfate soil",
        "target": 5,
        "queries": [
            'Australia AND "acid sulfate soil" AND (metagenomic OR metagenome)',
            'Oceania AND acidic soil AND (metagenomic OR metagenome)',
        ],
    },
    {
        "group": "Antarctica polar soil",
        "continent": "Antarctica",
        "acidic_soil_type": "Polar soil",
        "ecosystem": "Polar soil",
        "target": 5,
        "queries": [
            'Antarctica AND soil AND (metagenomic OR metagenome)',
            '"Antarctic soil" AND (metagenomic OR metagenome)',
        ],
    },
]


SOIL_POSITIVE = re.compile(
    r"\bsoil\b|rhizosphere|peat|peatland|bog|paddy|forest soil|red soil|acid sulfate|"
    r"permafrost|tundra|ferralsol|acrisol|oxisol|podzol|cropland|grassland",
    re.I,
)
NON_SOIL = re.compile(
    r"sediment|sludge|wastewater|groundwater|water column|lake water|marine|biofilm|"
    r"biofilter|bioreactor|rock metagenome|stromatolite|hot spring fluid|human|gut|feces",
    re.I,
)
ACIDIC_HINT = re.compile(
    r"acidic|acid soil|acid sulfate|red soil|paddy|peat|bog|permafrost|tundra|"
    r"ferralsol|acrisol|oxisol|podzol|Antarctic|polar",
    re.I,
)


def request_text(url: str, params: dict[str, str], sleep: float, retries: int = 4) -> str:
    query = urllib.parse.urlencode(params)
    last_error = None
    for attempt in range(1, retries + 1):
        time.sleep(sleep * attempt)
        try:
            req = urllib.request.Request(
                f"{url}?{query}",
                headers={"User-Agent": "ca-solibacter-pmta-global-abundance/0.1"},
            )
            with urllib.request.urlopen(req, timeout=180) as resp:
                return resp.read().decode("utf-8", errors="replace")
        except Exception as exc:  # noqa: BLE001
            last_error = exc
            if attempt == retries:
                raise
    raise RuntimeError(f"request failed: {last_error}")


def esearch(term: str, retmax: int, sleep: float) -> list[str]:
    text = request_text(
        ESEARCH,
        {"db": "sra", "term": term, "retmode": "json", "retmax": str(retmax)},
        sleep,
    )
    import json

    data = json.loads(text)
    return data.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str], sleep: float) -> list[dict[str, str]]:
    if not ids:
        return []
    text = request_text(
        ESUMMARY,
        {"db": "sra", "id": ",".join(ids), "retmode": "json"},
        sleep,
    )
    import json

    data = json.loads(text)
    rows = []
    for uid in data.get("result", {}).get("uids", []):
        item = data["result"][uid]
        run = match_one(item.get("runs", ""), r'Run acc="([^"]+)"')
        bioproject = match_one(item.get("expxml", ""), r"Bioproject>([^<]+)<")
        biosample = match_one(item.get("expxml", ""), r"Biosample>([^<]+)<")
        rows.append(
            {
                "sra_uid": uid,
                "run": run,
                "bioproject": bioproject,
                "biosample": biosample,
                "title": item.get("title", ""),
                "summary": item.get("summary", ""),
            }
        )
    return rows


def fetch_biosamples(accessions: list[str], sleep: float, chunk_size: int = 25) -> dict[str, dict[str, str]]:
    accessions = [x for x in dict.fromkeys(accessions) if x]
    out = {}
    for i in range(0, len(accessions), chunk_size):
        chunk = accessions[i : i + chunk_size]
        text = request_text(
            EFETCH,
            {"db": "biosample", "id": ",".join(chunk), "retmode": "xml"},
            sleep,
        )
        root = ET.fromstring(text)
        for bs in root.findall(".//BioSample"):
            acc = bs.attrib.get("accession", "")
            attrs = {}
            for attr in bs.findall(".//Attribute"):
                name = attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or ""
                if name:
                    attrs[name] = (attr.text or "").strip()
            title = text_or_empty(bs.find("./Description/Title"))
            organism = ""
            org = bs.find("./Description/Organism")
            if org is not None:
                organism = org.attrib.get("taxonomy_name") or text_or_empty(org.find("./OrganismName"))
            attrs["biosample_title"] = title
            attrs["biosample_organism"] = organism
            attrs["biosample_accession"] = acc
            out[acc] = attrs
    return out


def text_or_empty(elem) -> str:
    return "" if elem is None or elem.text is None else elem.text.strip()


def match_one(text: str, pattern: str) -> str:
    m = re.search(pattern, text or "")
    return m.group(1) if m else ""


def parse_lat_lon(text: str) -> tuple[str, str]:
    m = re.search(r"([0-9.]+)\s*([NS])\s+([0-9.]+)\s*([EW])", text or "")
    if not m:
        return "", ""
    lat = float(m.group(1))
    lon = float(m.group(3))
    if m.group(2) == "S":
        lat = -lat
    if m.group(4) == "W":
        lon = -lon
    return str(lat), str(lon)


def source_text(row: dict[str, str], bs: dict[str, str]) -> str:
    fields = [
        row.get("title", ""),
        row.get("summary", ""),
        bs.get("biosample_title", ""),
        bs.get("isolation_source", ""),
        bs.get("env_broad_scale", ""),
        bs.get("env_local_scale", ""),
        bs.get("env_medium", ""),
        bs.get("geo_loc_name", ""),
        bs.get("lat_lon", ""),
        bs.get("soil_type", ""),
        bs.get("description", ""),
    ]
    return "; ".join(x for x in fields if x)


def location_matches_continent(location: str, continent: str) -> bool:
    if not location:
        return True
    loc = location.lower()
    continent_terms = {
        "Asia": [
            "china",
            "taiwan",
            "japan",
            "korea",
            "india",
            "thailand",
            "vietnam",
            "indonesia",
            "malaysia",
            "philippines",
        ],
        "Europe": [
            "czech",
            "germany",
            "france",
            "spain",
            "italy",
            "sweden",
            "finland",
            "norway",
            "uk",
            "united kingdom",
            "poland",
            "netherlands",
            "europe",
        ],
        "North America": [
            "usa",
            "united states",
            "canada",
            "alaska",
            "minnesota",
            "wisconsin",
            "florida",
            "california",
            "north america",
        ],
        "South America": [
            "brazil",
            "peru",
            "colombia",
            "ecuador",
            "amazon",
            "argentina",
            "chile",
            "south america",
        ],
        "Africa": [
            "africa",
            "south africa",
            "kenya",
            "ghana",
            "nigeria",
            "tanzania",
            "ethiopia",
            "uganda",
            "limpopo",
        ],
        "Oceania": ["australia", "new zealand", "oceania"],
        "Antarctica": ["antarctica", "antarctic"],
    }
    terms = continent_terms.get(continent, [])
    return any(term in loc for term in terms)


def classify(row: dict[str, str], bs: dict[str, str], group: dict[str, object]) -> dict[str, str]:
    st = source_text(row, bs)
    is_soil = bool(SOIL_POSITIVE.search(st)) and not bool(NON_SOIL.search(st))
    acidic = bool(ACIDIC_HINT.search(st + " " + str(group["acidic_soil_type"])))
    lat, lon = parse_lat_lon(bs.get("lat_lon", ""))
    coord_status = "confirmed_lat_lon" if lat and lon else (
        "needs_geocoding" if bs.get("geo_loc_name") else "missing_location"
    )
    continent_ok = location_matches_continent(bs.get("geo_loc_name", ""), str(group["continent"]))
    if not is_soil:
        flag = "exclude_non_soil_or_unclear_source"
    elif not continent_ok:
        flag = "exclude_location_outside_continent"
    elif coord_status == "missing_location":
        flag = "exclude_missing_location"
    elif not acidic:
        flag = "candidate_pH_needs_confirmation"
    else:
        flag = "candidate"
    ph = bs.get("pH") or bs.get("ph") or bs.get("soil_pH") or bs.get("soil ph") or ""
    ph_status = "metadata_pH_available" if ph else "needs_metadata_or_OpenLandMap_confirmation"
    location = bs.get("geo_loc_name", "")
    return {
        "run": row["run"],
        "bioproject": row["bioproject"],
        "biosample": row["biosample"],
        "continent": str(group["continent"]),
        "acidic_soil_type": str(group["acidic_soil_type"]),
        "ecosystem": str(group["ecosystem"]),
        "location": location,
        "latitude": lat,
        "longitude": lon,
        "pH": ph,
        "pH_status": ph_status,
        "selection_reason": f"{group['continent']} representative; {group['acidic_soil_type']}",
        "project_dedup_group": row["bioproject"],
        "site_dedup_group": location or row["biosample"],
        "selection_rank": "",
        "coordinate_status": coord_status,
        "metadata_quality_flag": flag,
        "include_in_global_balanced_pilot": "no",
        "source_type": "external_acidic_soil_project_reads",
        "soil_only": "yes" if is_soil else "no",
        "acidic_soil": "yes" if acidic else "unknown",
        "soil_source_evidence": st,
    }


def select_group(candidates: list[dict[str, str]], target: int, max_per_project: int) -> list[dict[str, str]]:
    def score(row):
        coord_score = 0 if row["coordinate_status"] == "confirmed_lat_lon" else 1
        meta_score = 0 if row["metadata_quality_flag"] == "candidate" else 1
        return (coord_score, meta_score, row["bioproject"], row["biosample"], row["run"])

    selected = []
    project_counts = defaultdict(int)
    site_counts = defaultdict(int)
    used_biosamples = set()
    for row in sorted(candidates, key=score):
        if row["metadata_quality_flag"].startswith("exclude"):
            continue
        if row["biosample"] in used_biosamples:
            continue
        if project_counts[row["bioproject"]] >= max_per_project:
            continue
        if site_counts[row["site_dedup_group"]] >= 2:
            continue
        selected.append(row)
        used_biosamples.add(row["biosample"])
        project_counts[row["bioproject"]] += 1
        site_counts[row["site_dedup_group"]] += 1
        if len(selected) >= target:
            break
    for i, row in enumerate(selected, 1):
        row["selection_rank"] = str(i)
        row["include_in_global_balanced_pilot"] = "yes"
    return selected


def write_tsv(path: str, rows: list[dict[str, str]], fields: list[str]):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retmax", type=int, default=120)
    parser.add_argument("--sleep", type=float, default=0.45)
    parser.add_argument("--max-per-project", type=int, default=2)
    parser.add_argument("--out", default="results/global_acidic_soil_balanced_candidates_v2.tsv")
    parser.add_argument("--excluded-out", default="results/global_acidic_soil_balanced_excluded_v2.tsv")
    parser.add_argument("--accessions-out", default="config/global_acidic_soil_balanced_pilot_sra_accessions_v2.txt")
    args = parser.parse_args()

    all_candidates = []
    selected = []
    for group in QUERY_GROUPS:
        print(f"[group] {group['group']}")
        ids = []
        for query in group["queries"]:
            ids.extend(esearch(query, args.retmax, args.sleep))
        ids = list(dict.fromkeys(ids))
        summaries = []
        for i in range(0, len(ids), 100):
            summaries.extend(esummary(ids[i : i + 100], args.sleep))
        summaries = [row for row in summaries if row.get("run") and row.get("biosample") and row.get("bioproject")]
        biosamples = fetch_biosamples([row["biosample"] for row in summaries], args.sleep)
        classified = [classify(row, biosamples.get(row["biosample"], {}), group) for row in summaries]
        all_candidates.extend(classified)
        selected.extend(select_group(classified, int(group["target"]), args.max_per_project))

    seen_runs = set()
    unique_selected = []
    for row in selected:
        if row["run"] in seen_runs:
            continue
        seen_runs.add(row["run"])
        unique_selected.append(row)

    selected_runs = {row["run"] for row in unique_selected}
    excluded = []
    for row in all_candidates:
        if row["run"] in selected_runs:
            continue
        if row["metadata_quality_flag"].startswith("exclude"):
            excluded.append(row)
        elif row["include_in_global_balanced_pilot"] != "yes":
            row = dict(row)
            row["metadata_quality_flag"] = f"not_selected_after_project_dedup:{row['metadata_quality_flag']}"
            excluded.append(row)

    fields = [
        "run",
        "bioproject",
        "biosample",
        "continent",
        "acidic_soil_type",
        "ecosystem",
        "location",
        "latitude",
        "longitude",
        "pH",
        "pH_status",
        "selection_reason",
        "project_dedup_group",
        "site_dedup_group",
        "selection_rank",
        "coordinate_status",
        "metadata_quality_flag",
        "include_in_global_balanced_pilot",
        "source_type",
        "soil_only",
        "acidic_soil",
        "soil_source_evidence",
    ]
    write_tsv(args.out, unique_selected, fields)
    write_tsv(args.excluded_out, excluded, fields)
    Path(args.accessions_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.accessions_out).write_text(
        "\n".join(row["run"] for row in unique_selected) + ("\n" if unique_selected else ""),
        encoding="ascii",
    )
    print(f"Wrote {len(unique_selected)} selected rows to {args.out}")
    print(f"Wrote {len(excluded)} excluded rows to {args.excluded_out}")
    print(f"Wrote accession list to {args.accessions_out}")


if __name__ == "__main__":
    main()
