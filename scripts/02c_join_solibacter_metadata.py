from pathlib import Path
import csv
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = ROOT / "results" / "candidate_solibacter_manifest.tsv"
OUT = ROOT / "results" / "solibacter_sample_manifest.tsv"

REL_NS = [
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "http://purl.oclc.org/ooxml/officeDocument/relationships",
]
PKG_NS = {"pkg": "http://schemas.openxmlformats.org/package/2006/relationships"}


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


def descendants_by_name(node, name):
    return [child for child in node.iter() if local_name(child.tag) == name]


def children_by_name(node, name):
    return [child for child in list(node) if local_name(child.tag) == name]


def column_index(cell_ref):
    letters = re.sub(r"[^A-Z]", "", cell_ref.upper())
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch) - ord("A") + 1)
    return idx - 1


def shared_strings(zf):
    if "xl/sharedStrings.xml" not in zf.namelist():
        return []
    root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
    return ["".join(t.text or "" for t in descendants_by_name(si, "t")) for si in descendants_by_name(root, "si")]


def workbook_sheets(zf):
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels.findall("pkg:Relationship", PKG_NS)}
    for sheet in descendants_by_name(wb, "sheet"):
        rid = None
        for ns in REL_NS:
            rid = sheet.attrib.get(f"{{{ns}}}id")
            if rid:
                break
        if rid in rel_map:
            yield sheet.attrib.get("name", "sheet"), "xl/" + rel_map[rid].lstrip("/")


def cell_value(cell, strings):
    ctype = cell.attrib.get("t")
    values = children_by_name(cell, "v")
    if not values:
        return ""
    raw = values[0].text or ""
    if ctype == "s" and raw.isdigit():
        return strings[int(raw)]
    return raw


def read_xlsx_first_sheet(path):
    with zipfile.ZipFile(path) as zf:
        strings = shared_strings(zf)
        _, sheet_path = next(workbook_sheets(zf))
        root = ET.fromstring(zf.read(sheet_path))
        rows = []
        for row in descendants_by_name(root, "row"):
            values = {}
            max_idx = -1
            for cell in children_by_name(row, "c"):
                idx = column_index(cell.attrib.get("r", "A1"))
                values[idx] = cell_value(cell, strings)
                max_idx = max(max_idx, idx)
            rows.append([values.get(i, "") for i in range(max_idx + 1)])
        header_idx = 1 if rows and rows[0] and str(rows[0][0]).startswith("Supplementary Data") else 0
        header = rows[header_idx]
        out = []
        for row in rows[header_idx + 1 :]:
            padded = row + [""] * (len(header) - len(row))
            out.append(dict(zip(header, padded)))
        return out


def sample_prefix(mag_id):
    return re.sub(r"\.[0-9]+$", "", mag_id)


if not MANIFEST.exists():
    print("Run 02b_make_solibacter_manifest.py first.", file=sys.stderr)
    sys.exit(1)

metadata_rows = read_xlsx_first_sheet(RAW / "Supplementary Data 1.xlsx")
metadata_by_user_genome = {row.get("user_genome", ""): row for row in metadata_rows}

with MANIFEST.open(encoding="utf-8", newline="") as fh:
    candidates = list(csv.DictReader(fh, delimiter="\t"))

fields = [
    "mag_id",
    "sample_prefix",
    "sample",
    "project",
    "source",
    "nation",
    "continent",
    "ecosystem",
    "paired_seqs",
    "latitude",
    "longitude",
    "bioproject",
    "biosample",
    "soil_pH",
    "metadata_status",
]

with OUT.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for cand in candidates:
        prefix = sample_prefix(cand["mag_id"])
        meta = metadata_by_user_genome.get(prefix)
        if meta is None:
            writer.writerow({"mag_id": cand["mag_id"], "sample_prefix": prefix, "metadata_status": "missing"})
            continue
        writer.writerow(
            {
                "mag_id": cand["mag_id"],
                "sample_prefix": prefix,
                "sample": meta.get("Sample", ""),
                "project": meta.get("Project", ""),
                "source": meta.get("Source", ""),
                "nation": meta.get("Nation", ""),
                "continent": meta.get("Continent", ""),
                "ecosystem": meta.get("Ecosystem", ""),
                "paired_seqs": meta.get("Paired_seqs", ""),
                "latitude": meta.get("Lat", ""),
                "longitude": meta.get("Lon", ""),
                "bioproject": meta.get("Bioproject", ""),
                "biosample": meta.get("Biosample", ""),
                "soil_pH": "",
                "metadata_status": "matched_no_pH_in_supplementary_data_1",
            }
        )

print(f"Wrote metadata join table to {OUT}")
