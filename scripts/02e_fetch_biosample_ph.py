from pathlib import Path
import csv
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "results" / "solibacter_sample_manifest.tsv"
OUT = ROOT / "results" / "biosample_ph_recovered.tsv"

PH_KEYS = ("ph", "soil ph", "soil_ph", "env ph", "env_ph")


def fetch_biosample_xml(accession):
    term = urllib.parse.quote(accession)
    esearch = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=biosample&term={term}&retmode=xml"
    with urllib.request.urlopen(esearch, timeout=30) as fh:
        root = ET.fromstring(fh.read())
    ids = [node.text for node in root.findall(".//Id") if node.text]
    if not ids:
        return None
    efetch = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id={ids[0]}&retmode=xml"
    with urllib.request.urlopen(efetch, timeout=30) as fh:
        return fh.read()


def parse_attributes(xml_bytes):
    attrs = {}
    if not xml_bytes:
        return attrs
    root = ET.fromstring(xml_bytes)
    for attr in root.findall(".//Attribute"):
        name = attr.attrib.get("attribute_name") or attr.attrib.get("harmonized_name") or ""
        value = (attr.text or "").strip()
        if name or value:
            attrs[name] = value
    return attrs


def find_ph(attrs):
    for key, value in attrs.items():
        norm = key.strip().lower()
        if norm in PH_KEYS or norm.replace(" ", "_") in PH_KEYS:
            return value, key
        if "ph" in norm:
            match = re.search(r"([0-9]+(?:\.[0-9]+)?)", value)
            if match:
                ph = float(match.group(1))
                if 2.0 <= ph <= 12.0:
                    return match.group(1), key
    return "", ""


with META.open(encoding="utf-8", newline="") as fh:
    rows = list(csv.DictReader(fh, delimiter="\t"))

biosamples = sorted({row["biosample"] for row in rows if row.get("biosample")})

with OUT.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(
        fh,
        fieldnames=["biosample", "soil_pH", "pH_source_key", "status", "attribute_count"],
        delimiter="\t",
    )
    writer.writeheader()
    for i, biosample in enumerate(biosamples, start=1):
        try:
            xml_bytes = fetch_biosample_xml(biosample)
            attrs = parse_attributes(xml_bytes)
            ph, source = find_ph(attrs)
            writer.writerow(
                {
                    "biosample": biosample,
                    "soil_pH": ph,
                    "pH_source_key": source,
                    "status": "pH_found" if ph else "pH_not_found",
                    "attribute_count": len(attrs),
                }
            )
        except Exception as exc:
            writer.writerow(
                {
                    "biosample": biosample,
                    "soil_pH": "",
                    "pH_source_key": "",
                    "status": f"fetch_error: {exc}",
                    "attribute_count": 0,
                }
            )
        if i % 3 == 0:
            time.sleep(0.4)

print(f"Wrote BioSample pH recovery table to {OUT}")
