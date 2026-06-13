from pathlib import Path
import csv
import json
import math
import subprocess
import urllib.request
import shutil
import tempfile

ROOT = Path(__file__).resolve().parents[1]
IN_HITS = ROOT / "results" / "metagraph" / "metagraph_hits_limit500.tsv"
IN_SUMMARY = ROOT / "results" / "metagraph" / "metagraph_result_summary_limit500.tsv"
OUT_DIR = ROOT / "results" / "metagraph" / "figures"
WORLD = ROOT / "data" / "raw" / "world.geojson"
OUT_DIR.mkdir(parents=True, exist_ok=True)
WORLD.parent.mkdir(parents=True, exist_ok=True)

WIDTH = 1200
HEIGHT = 650
MAP_LEFT = 70
MAP_TOP = 55
MAP_WIDTH = 1060
MAP_HEIGHT = 470


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fnum(value):
    try:
        return float(value)
    except Exception:
        return None


def esc(text):
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def project(lon, lat):
    x = MAP_LEFT + (lon + 180.0) / 360.0 * MAP_WIDTH
    y = MAP_TOP + (90.0 - lat) / 180.0 * MAP_HEIGHT
    return x, y


def download_world():
    if WORLD.exists() and WORLD.stat().st_size > 1000:
        return
    url = "https://raw.githubusercontent.com/holtzy/D3-graph-gallery/master/DATA/world.geojson"
    with urllib.request.urlopen(url, timeout=60) as resp:
        WORLD.write_bytes(resp.read())


def polygon_path(coords):
    parts = []
    prev_x = None
    for ring in coords:
        if not ring:
            continue
        cmds = []
        ring_prev_x = None
        broken = False
        for lon, lat, *_ in ring:
            x, y = project(float(lon), float(lat))
            if ring_prev_x is not None and abs(x - ring_prev_x) > MAP_WIDTH * 0.55:
                broken = True
                break
            cmds.append((x, y))
            ring_prev_x = x
        if broken or not cmds:
            continue
        first = cmds[0]
        parts.append(f"M{first[0]:.1f},{first[1]:.1f}")
        for x, y in cmds[1:]:
            parts.append(f"L{x:.1f},{y:.1f}")
        parts.append("Z")
        prev_x = cmds[-1][0]
    return " ".join(parts)


def world_paths():
    try:
        download_world()
        data = json.loads(WORLD.read_text(encoding="utf-8"))
    except Exception:
        return ""
    paths = []
    for feature in data.get("features", []):
        geom = feature.get("geometry") or {}
        gtype = geom.get("type")
        coords = geom.get("coordinates") or []
        if gtype == "Polygon":
            d = polygon_path(coords)
            if d:
                paths.append(f'<path d="{d}" class="land"/>')
        elif gtype == "MultiPolygon":
            for poly in coords:
                d = polygon_path(poly)
                if d:
                    paths.append(f'<path d="{d}" class="land"/>')
    return "\n".join(paths)


def color(score):
    # Blue to orange-red gradient.
    s = max(0.0, min(1.0, score))
    stops = [
        (0.0, (49, 104, 142)),
        (0.35, (53, 161, 151)),
        (0.7, (242, 183, 5)),
        (1.0, (211, 64, 53)),
    ]
    for i in range(len(stops) - 1):
        a, ca = stops[i]
        b, cb = stops[i + 1]
        if a <= s <= b:
            t = (s - a) / (b - a) if b > a else 0
            rgb = tuple(round(ca[j] + (cb[j] - ca[j]) * t) for j in range(3))
            return f"rgb({rgb[0]},{rgb[1]},{rgb[2]})"
    return "rgb(211,64,53)"


def graticule():
    lines = []
    for lon in range(-180, 181, 60):
        x1, y1 = project(lon, -60)
        x2, y2 = project(lon, 85)
        lines.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="grid"/>')
        lines.append(f'<text x="{x1:.1f}" y="{MAP_TOP + MAP_HEIGHT + 24}" class="axis">{lon}</text>')
    for lat in range(-60, 91, 30):
        x1, y1 = project(-180, lat)
        x2, y2 = project(180, lat)
        lines.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="grid"/>')
        lines.append(f'<text x="{MAP_LEFT - 18}" y="{y1 + 4:.1f}" class="axis" text-anchor="end">{lat}</text>')
    return "\n".join(lines)


def map_svg(rows, database, title, out_stem):
    rows = [r for r in rows if r["database"] == database]
    geo = []
    for r in rows:
        lat = fnum(r.get("latitude"))
        lon = fnum(r.get("longitude"))
        score = fnum(r.get("normalized_score"))
        kmer = fnum(r.get("kmer_count"))
        if lat is None or lon is None or score is None or kmer is None:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        geo.append((lat, lon, score, kmer, r))
    max_kmer = max([g[3] for g in geo], default=1)
    circles = []
    for lat, lon, score, kmer, r in geo:
        x, y = project(lon, lat)
        radius = 2.5 + 9.5 * math.sqrt(kmer / max_kmer)
        circles.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.2f}" '
            f'fill="{color(score)}" fill-opacity="0.68" stroke="#253044" stroke-opacity="0.35" stroke-width="0.45">'
            f'<title>{esc(r.get("sequence_id",""))} | {esc(r.get("location",""))} | normalized score={score:.3f} | k-mer count={int(kmer)}</title>'
            f"</circle>"
        )
    legend = []
    for i, val in enumerate([0.2, 0.4, 0.6, 0.8, 1.0]):
        x = 835 + i * 48
        legend.append(f'<rect x="{x}" y="578" width="42" height="12" fill="{color(val)}"/>')
        legend.append(f'<text x="{x+21}" y="607" class="legend" text-anchor="middle">{val:.1f}</text>')
    size_legend = []
    for i, frac in enumerate([0.25, 0.55, 1.0]):
        r = 2.5 + 9.5 * math.sqrt(frac)
        x = 120 + i * 70
        size_legend.append(f'<circle cx="{x}" cy="586" r="{r:.2f}" fill="#7699c9" fill-opacity="0.55" stroke="#253044" stroke-width="0.5"/>')
        size_legend.append(f'<text x="{x}" y="622" class="legend" text-anchor="middle">{int(max_kmer*frac)}</text>')
    subtitle = (
        f"{len(geo)} geolocated MetaGraph preview hits; point colour = normalized k-mer score; "
        "point size = k-mer count. Search signal, not quantitative abundance."
    )
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">
<style>
  .bg {{ fill: #f8fafc; }}
  .land {{ fill: #e6e2d8; stroke: #b8b2a7; stroke-width: 0.45; }}
  .grid {{ stroke: #cbd5e1; stroke-width: 0.55; stroke-dasharray: 3 5; }}
  .title {{ font-family: Arial, sans-serif; font-size: 24px; font-weight: 700; fill: #172033; }}
  .subtitle {{ font-family: Arial, sans-serif; font-size: 13px; fill: #475569; }}
  .axis {{ font-family: Arial, sans-serif; font-size: 10px; fill: #64748b; }}
  .legend {{ font-family: Arial, sans-serif; font-size: 11px; fill: #334155; }}
  .note {{ font-family: Arial, sans-serif; font-size: 12px; fill: #475569; }}
</style>
<rect width="100%" height="100%" class="bg"/>
<text x="60" y="32" class="title">{esc(title)}</text>
<text x="60" y="52" class="subtitle">{esc(subtitle)}</text>
<rect x="{MAP_LEFT}" y="{MAP_TOP}" width="{MAP_WIDTH}" height="{MAP_HEIGHT}" fill="#eef6fb" stroke="#94a3b8" stroke-width="0.8"/>
{graticule()}
{world_paths()}
{"".join(circles)}
<text x="90" y="560" class="legend">k-mer count</text>
{"".join(size_legend)}
<text x="835" y="560" class="legend">normalized k-mer score</text>
{"".join(legend)}
<text x="60" y="640" class="note">Exact Ellin6076 pmtA nucleotide and translated PmtA searches returned no hits in tested MetaGraph databases.</text>
</svg>"""
    save_svg_html_pdf(svg, out_stem)
    return len(geo)


def bar_svg(rows, out_stem):
    grouped = {}
    for database in ["sra-microbe", "refseq33m"]:
        sub = [r for r in rows if r["database"] == database and r.get("location")]
        counts = {}
        for r in sub:
            loc = r["location"].strip() or "Unknown"
            counts[loc] = counts.get(loc, 0) + 1
        grouped[database] = sorted(counts.items(), key=lambda x: x[1], reverse=True)[:10]
    bar_w = 36
    panel_w = 500
    max_count = max([v for items in grouped.values() for _, v in items], default=1)
    panels = []
    for p, database in enumerate(["sra-microbe", "refseq33m"]):
        x0 = 80 + p * 560
        y0 = 540
        panels.append(f'<text x="{x0}" y="82" class="panel">{database}</text>')
        panels.append(f'<line x1="{x0}" y1="{y0}" x2="{x0+panel_w}" y2="{y0}" class="axisline"/>')
        panels.append(f'<line x1="{x0}" y1="115" x2="{x0}" y2="{y0}" class="axisline"/>')
        for i, (loc, count) in enumerate(grouped[database]):
            x = x0 + 20 + i * 45
            h = 390 * count / max_count
            y = y0 - h
            panels.append(f'<rect x="{x}" y="{y:.1f}" width="{bar_w}" height="{h:.1f}" fill="#4f8fbf"/>')
            panels.append(f'<text x="{x+bar_w/2}" y="{y-6:.1f}" class="count" text-anchor="middle">{count}</text>')
            label = loc if len(loc) <= 16 else loc[:15] + "..."
            panels.append(f'<text x="{x+bar_w/2}" y="{y0+14}" class="xlabel" transform="rotate(55 {x+bar_w/2},{y0+14})">{esc(label)}</text>')
    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">
<style>
  .bg {{ fill: #ffffff; }}
  .title {{ font-family: Arial, sans-serif; font-size: 24px; font-weight: 700; fill: #172033; }}
  .subtitle {{ font-family: Arial, sans-serif; font-size: 13px; fill: #475569; }}
  .panel {{ font-family: Arial, sans-serif; font-size: 16px; font-weight: 700; fill: #172033; }}
  .axisline {{ stroke: #334155; stroke-width: 1; }}
  .xlabel {{ font-family: Arial, sans-serif; font-size: 10px; fill: #334155; text-anchor: start; }}
  .count {{ font-family: Arial, sans-serif; font-size: 10px; fill: #334155; }}
</style>
<rect width="100%" height="100%" class="bg"/>
<text x="60" y="38" class="title">Location-level Ellin6076-related 16S detection intensity</text>
<text x="60" y="60" class="subtitle">Top 10 location labels by MetaGraph preview hit count. Counts represent search-detection records, not quantitative abundance.</text>
{"".join(panels)}
</svg>"""
    save_svg_html_pdf(svg, out_stem)


def save_svg_html_pdf(svg, out_stem):
    svg_path = OUT_DIR / f"{out_stem}.svg"
    html_path = OUT_DIR / f"{out_stem}.html"
    pdf_path = OUT_DIR / f"{out_stem}.pdf"
    svg_path.write_text(svg, encoding="utf-8")
    html = f"""<!doctype html><html><head><meta charset="utf-8"><style>
@page {{ size: {WIDTH}px {HEIGHT}px; margin: 0; }}
html, body {{ margin: 0; padding: 0; width: {WIDTH}px; height: {HEIGHT}px; }}
</style></head><body>{svg}</body></html>"""
    html_path.write_text(html, encoding="utf-8")
    chrome_paths = [
        Path(r"C:\Program Files\Google\Chrome\Application\chrome.exe"),
        Path(r"C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe"),
    ]
    exe = next((p for p in chrome_paths if p.exists()), None)
    if exe:
        try:
            tmp_dir = Path(tempfile.gettempdir()) / "metagraph_pdf_export"
            tmp_dir.mkdir(parents=True, exist_ok=True)
            tmp_html = tmp_dir / f"{out_stem}.html"
            tmp_pdf = tmp_dir / f"{out_stem}.pdf"
            tmp_html.write_text(html, encoding="utf-8")
            if tmp_pdf.exists():
                tmp_pdf.unlink()
            subprocess.run(
                [
                    str(exe),
                    "--headless=new",
                    "--no-sandbox",
                    "--disable-dev-shm-usage",
                    "--disable-gpu",
                    f"--print-to-pdf={tmp_pdf}",
                    f"file:///{tmp_html.resolve().as_posix()}",
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=60,
            )
            if tmp_pdf.exists():
                shutil.copy2(tmp_pdf, pdf_path)
            err_path = OUT_DIR / f"{out_stem}_pdf_error.txt"
            if err_path.exists():
                err_path.unlink()
        except Exception as exc:
            (OUT_DIR / f"{out_stem}_pdf_error.txt").write_text(str(exc), encoding="utf-8")


def write_summary(rows):
    valid = []
    for r in rows:
        lat = fnum(r.get("latitude"))
        lon = fnum(r.get("longitude"))
        score = fnum(r.get("normalized_score"))
        kmer = fnum(r.get("kmer_count"))
        if lat is None or lon is None or score is None or kmer is None:
            continue
        if not (-90 <= lat <= 90 and -180 <= lon <= 180):
            continue
        valid.append((r, score, kmer))
    out = OUT_DIR / "Ellin6076_metagraph_distribution_summary.tsv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["database", "total_hits", "geolocated_hits", "max_normalized_score", "mean_normalized_score", "max_kmer_count", "note"])
        for database in ["sra-microbe", "refseq33m", "atb"]:
            total = sum(1 for r in rows if r["database"] == database)
            sub = [(score, kmer) for r, score, kmer in valid if r["database"] == database]
            scores = [s for s, _ in sub]
            kmers = [k for _, k in sub]
            writer.writerow([
                database,
                total,
                len(sub),
                f"{max(scores):.6f}" if scores else "",
                f"{sum(scores)/len(scores):.6f}" if scores else "",
                int(max(kmers)) if kmers else "",
                "MetaGraph 16S search signal; not quantitative abundance",
            ])


def main():
    rows = [r for r in read_tsv(IN_HITS) if r["query"] == "Ellin6076_16S_partial"]
    write_summary(rows)
    map_svg(rows, "sra-microbe", "Global distribution of Ellin6076-related 16S search signals (SRA microbe)", "Ellin6076_16S_global_distribution_sra_microbe")
    map_svg(rows, "refseq33m", "Global distribution of Ellin6076-related 16S search signals (RefSeq)", "Ellin6076_16S_global_distribution_refseq33m")
    bar_svg(rows, "Ellin6076_16S_detection_intensity_by_location")
    print(f"Wrote figures to {OUT_DIR}")


if __name__ == "__main__":
    main()
