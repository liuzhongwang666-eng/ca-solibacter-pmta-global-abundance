#!/usr/bin/env python3
import argparse
import csv
import gzip
from pathlib import Path


def read_fasta_lengths(path):
    lengths = {}
    name = None
    seqs = []
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = len("".join(seqs).replace("*", ""))
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
    if name is not None:
        lengths[name] = len("".join(seqs).replace("*", ""))
    return lengths


def count_fastq_reads(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    n_lines = 0
    with opener(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for _ in fh:
            n_lines += 1
    return n_lines // 4


def best_filtered_hits(path, min_identity, min_alignment_aa):
    best = {}
    if not Path(path).exists() or Path(path).stat().st_size == 0:
        return best
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            length = int(parts[3])
            evalue = float(parts[10])
            bitscore = float(parts[11])
            if pident < min_identity or length < min_alignment_aa:
                continue
            current = best.get(qseqid)
            if current is None or bitscore > current["bitscore"] or (
                bitscore == current["bitscore"] and evalue < current["evalue"]
            ):
                best[qseqid] = {
                    "sseqid": sseqid,
                    "pident": pident,
                    "alignment_length_aa": length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                }
    return best


parser = argparse.ArgumentParser(
    description="Summarize DIAMOND blastx hits to functional PmtA marker abundance."
)
parser.add_argument("--blastx-dir", required=True)
parser.add_argument("--reads-dir", required=True)
parser.add_argument("--samples", required=True)
parser.add_argument("--markers", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--min-identity", type=float, default=35.0)
parser.add_argument("--min-alignment-aa", type=int, default=30)
args = parser.parse_args()

blastx_dir = Path(args.blastx_dir)
reads_dir = Path(args.reads_dir)
samples = [line.strip() for line in Path(args.samples).read_text(encoding="utf-8").splitlines() if line.strip()]
marker_lengths = read_fasta_lengths(args.markers)
total_marker_aa_kb = sum(marker_lengths.values()) / 1000 if marker_lengths else 0

rows = []
for sample in samples:
    r1 = reads_dir / f"{sample}_1.fastq.gz"
    r2 = reads_dir / f"{sample}_2.fastq.gz"
    total_reads = count_fastq_reads(r1) + count_fastq_reads(r2)
    sample_best = {}
    marker_counts = {marker: 0 for marker in marker_lengths}

    for mate in (1, 2):
        hit_file = blastx_dir / f"{sample}_{mate}.pmta_marker.blastx.tsv"
        hits = best_filtered_hits(hit_file, args.min_identity, args.min_alignment_aa)
        for read_id, hit in hits.items():
            key = f"{mate}:{read_id}"
            sample_best[key] = hit

    for hit in sample_best.values():
        marker = hit["sseqid"]
        marker_counts[marker] = marker_counts.get(marker, 0) + 1

    hit_reads = sum(marker_counts.values())
    reads_per_million = hit_reads / (total_reads / 1_000_000) if total_reads else 0
    rpkm_like = (
        hit_reads / (total_marker_aa_kb * (total_reads / 1_000_000))
        if total_marker_aa_kb and total_reads
        else 0
    )
    top_marker = ""
    top_marker_hits = 0
    if marker_counts:
        top_marker, top_marker_hits = max(marker_counts.items(), key=lambda item: item[1])

    rows.append(
        {
            "sample_id": sample,
            "total_reads_searched": total_reads,
            "pmta_marker_hit_reads": hit_reads,
            "pmta_marker_reads_per_million": reads_per_million,
            "pmta_marker_rpkm_like": rpkm_like,
            "n_markers_in_database": len(marker_lengths),
            "n_markers_with_hits": sum(1 for v in marker_counts.values() if v > 0),
            "top_marker": top_marker,
            "top_marker_hits": top_marker_hits,
            "min_identity": args.min_identity,
            "min_alignment_aa": args.min_alignment_aa,
        }
    )

out = Path(args.out)
out.parent.mkdir(parents=True, exist_ok=True)
fieldnames = [
    "sample_id",
    "total_reads_searched",
    "pmta_marker_hit_reads",
    "pmta_marker_reads_per_million",
    "pmta_marker_rpkm_like",
    "n_markers_in_database",
    "n_markers_with_hits",
    "top_marker",
    "top_marker_hits",
    "min_identity",
    "min_alignment_aa",
]
with out.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} rows to {out}")
