from pathlib import Path
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "results" / "candidate_solibacter_mags.tsv"
PATTERNS = (
    re.compile(r"\bcandidatus\s+solibacter\b", re.I),
    re.compile(r"\bca\.\s*solibacter\b", re.I),
    re.compile(r"\bg__solibacter\b", re.I),
    re.compile(r"\bs__solibacter\b", re.I),
    re.compile(r"\bf__solibacteraceae\b", re.I),
)
REL_NS = [
    "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    "http://purl.oclc.org/ooxml/officeDocument/relationships",
]
PKG_NS = {"pkg": "http://schemas.openxmlformats.org/package/2006/relationships"}


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


def children_by_name(node, name):
    return [child for child in list(node) if local_name(child.tag) == name]


def descendants_by_name(node, name):
    return [child for child in node.iter() if local_name(child.tag) == name]


def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


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
    values = []
    for si in descendants_by_name(root, "si"):
        texts = [t.text or "" for t in descendants_by_name(si, "t")]
        values.append("".join(texts))
    return values


def workbook_sheets(zf):
    wb = ET.fromstring(zf.read("xl/workbook.xml"))
    rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
    rel_map = {
        rel.attrib["Id"]: rel.attrib["Target"]
        for rel in rels.findall("pkg:Relationship", PKG_NS)
    }
    sheets = []
    for sheet in descendants_by_name(wb, "sheet"):
        name = sheet.attrib.get("name", "sheet")
        rid = None
        for ns in REL_NS:
            rid = sheet.attrib.get(f"{{{ns}}}id")
            if rid:
                break
        target = rel_map.get(rid)
        if target:
            sheets.append((name, "xl/" + target.lstrip("/")))
    return sheets


def cell_value(cell, strings):
    ctype = cell.attrib.get("t")
    values = children_by_name(cell, "v")
    v = values[0] if values else None
    if ctype == "inlineStr":
        texts = [t.text or "" for t in descendants_by_name(cell, "t")]
        return "".join(texts)
    if v is None:
        return ""
    raw = v.text or ""
    if ctype == "s":
        try:
            return strings[int(raw)]
        except Exception:
            return raw
    return raw


def iter_rows(zf, sheet_path, strings):
    root = ET.fromstring(zf.read(sheet_path))
    for row in descendants_by_name(root, "row"):
        values = {}
        max_idx = -1
        for cell in children_by_name(row, "c"):
            idx = column_index(cell.attrib.get("r", "A1"))
            values[idx] = cell_value(cell, strings)
            max_idx = max(max_idx, idx)
        yield [values.get(i, "") for i in range(max_idx + 1)]


def clean_header(value, idx):
    value = str(value).strip()
    if not value:
        value = f"unnamed_{idx+1}"
    value = re.sub(r"[^0-9A-Za-z_.]+", "_", value).strip("_")
    if not value:
        value = f"unnamed_{idx+1}"
    return value


files = [
    RAW / "Supplementary Data 2.xlsx",
    RAW / "Supplementary Data 2 subgroup.xlsx",
]
existing = [p for p in files if p.exists()]
if not existing:
    fail("Put SMAG Supplementary Data 2.xlsx or Supplementary Data 2 subgroup.xlsx in data/raw/.")

out_rows = []
max_cols = 0
for path in existing:
    with zipfile.ZipFile(path) as zf:
        strings = shared_strings(zf)
        for sheet_name, sheet_path in workbook_sheets(zf):
            rows = list(iter_rows(zf, sheet_path, strings))
            if len(rows) < 2:
                continue
            header_idx = 1 if rows[0] and str(rows[0][0]).startswith("Supplementary Data") else 0
            header = [clean_header(v, i) for i, v in enumerate(rows[header_idx])]
            data_start = header_idx + 1
            for row_num, row in enumerate(rows[data_start:], start=data_start + 1):
                row_text = " | ".join(row).lower()
                if any(pattern.search(row_text) for pattern in PATTERNS):
                    max_cols = max(max_cols, len(header), len(row))
                    out_rows.append((path.name, sheet_name, row_num, header, row))

if not out_rows:
    OUT.write_text("source_file\tsource_sheet\trow_number\n", encoding="utf-8")
    print("No Candidatus Solibacter-like rows found. Check taxonomy terms.")
else:
    # Use the longest observed header, then fill any missing columns.
    best_header = max((header for _, _, _, header, _ in out_rows), key=len)
    best_header = best_header + [f"extra_col_{i+1}" for i in range(len(best_header), max_cols)]
    lines = ["source_file\tsource_sheet\trow_number\t" + "\t".join(best_header)]
    for source_file, sheet_name, row_num, header, row in out_rows:
        padded = row + [""] * (max_cols - len(row))
        safe = [str(v).replace("\t", " ").replace("\n", " ") for v in padded]
        lines.append(f"{source_file}\t{sheet_name}\t{row_num}\t" + "\t".join(safe))
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {len(out_rows)} candidate rows to {OUT}")
