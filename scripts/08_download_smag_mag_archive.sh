#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PART_DIR="${ROOT}/data/raw/smag_mag_parts"
OUT_ARCHIVE="${ROOT}/data/raw/mag.tar.gz"
MAG_DIR="${ROOT}/data/raw/mag"
ZENODO_RECORD="${ZENODO_RECORD:-8223844}"

mkdir -p "${PART_DIR}" "${MAG_DIR}"

echo "[1/4] Querying Zenodo record ${ZENODO_RECORD} for SMAG MAG archive parts..."
python3 - <<PY > "${PART_DIR}/mag_part_urls.tsv"
import json
import urllib.request

record = "${ZENODO_RECORD}"
url = f"https://zenodo.org/api/records/{record}"
with urllib.request.urlopen(url, timeout=60) as fh:
    data = json.load(fh)

rows = []
for item in data.get("files", []):
    key = item.get("key", "")
    links = item.get("links", {})
    download = links.get("self") or links.get("download")
    if key.startswith("mag.tar.gz."):
        rows.append((key, download))

if not rows:
    raise SystemExit("No mag.tar.gz.* files found in Zenodo record.")

for key, download in sorted(rows):
    print(f"{key}\t{download}")
PY

echo "[2/4] Downloading MAG archive parts..."
while IFS=$'\t' read -r name url; do
  dest="${PART_DIR}/${name}"
  if [[ -s "${dest}" ]]; then
    echo "Exists: ${dest}"
  else
    echo "Downloading: ${name}"
    if command -v wget >/dev/null 2>&1; then
      wget -c -O "${dest}" "${url}"
    else
      curl -L --retry 5 -C - -o "${dest}" "${url}"
    fi
  fi
done < "${PART_DIR}/mag_part_urls.tsv"

echo "[3/4] Combining MAG archive parts..."
cat "${PART_DIR}"/mag.tar.gz.* > "${OUT_ARCHIVE}"

echo "[4/4] Extracting MAG archive to ${MAG_DIR}..."
tar -xzf "${OUT_ARCHIVE}" -C "${MAG_DIR}"

echo "Done. MAG files are under: ${MAG_DIR}"
