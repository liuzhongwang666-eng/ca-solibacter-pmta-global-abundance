#!/usr/bin/env python3
"""Summarize acidic-soil extension Solibacter and PmtA marker outputs."""

import argparse
import csv
import gzip
from pathlib import Path


def read_table(path, delimiter="\t"):
    p = Path(path)
    if not p.exists():
        return []
    with p.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def find_col(fields, text):
    for field in fields:
        if text.lower() in field.lower():
            return field
    return None


def count_fastq_reads(path):
    p = Path(path)
    if not p.exists():
        return ""
    opener = gzip.open if str(p).endswith(".gz") else open
    n = 0
    with opener(p, "rt", encoding="utf-8", errors="ignore") as fh:
        for _ in fh:
            n += 1
    return n // 4


def marker_lengths(path):
    lengths = {}
    name = None
    seqs = []
    p = Path(path)
    if not p.exists():
        return lengths
    with p.open(encoding="utf-8") as fh:
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


def best_filtered_hits(path, min_identity, min_alignment_aa):
    best = {}
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return best
    with p.open(encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            qseqid, sseqid = parts[0], parts[1]
            pident, length = float(parts[2]), int(parts[3])
            evalue, bitscore = float(parts[10]), float(parts[11])
            if pident < min_identity or length < min_alignment_aa:
                continue
            current = best.get(qseqid)
            if current is None or bitscore > current["bitscore"] or (
                bitscore == current["bitscore"] and evalue < current["evalue"]
            ):
                best[qseqid] = {"sseqid": sseqid, "bitscore": bitscore, "evalue": evalue}
    return best


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", default="config/acidic_soil_pilot_sra_accessions.txt")
    parser.add_argument("--metadata", default="results/acidic_soil_sra_extension_candidates.tsv")
    parser.add_argument("--coverm-dir", default="results/coverm_solibacter_acidic_extension")
    parser.add_argument("--blastx-dir", default="results/pmta_marker_blastx_acidic_extension")
    parser.add_argument("--reads-dir", default="data/raw/reads")
    parser.add_argument("--markers", default="data/processed/pmta_markers/functional_pmta_marker_candidates.faa")
    parser.add_argument("--out", default="results/acidic_soil_solibacter_pmta_abundance_summary.csv")
    parser.add_argument("--min-identity", type=float, default=35.0)
    parser.add_argument("--min-alignment-aa", type=int, default=30)
    args = parser.parse_args()

    samples = [x.strip() for x in Path(args.samples).read_text(encoding="utf-8").splitlines() if x.strip()]
    metadata = {row.get("run", ""): row for row in read_table(args.metadata)}
    marker_len = marker_lengths(args.markers)
    marker_aa_kb = sum(marker_len.values()) / 1000 if marker_len else 0

    rows_out = []
    for sample in samples:
        coverm_path = Path(args.coverm_dir) / f"{sample}.coverm.tsv"
        r1 = Path(args.reads_dir) / f"{sample}_1.fastq.gz"
        r2 = Path(args.reads_dir) / f"{sample}_2.fastq.gz"
        b1 = Path(args.blastx_dir) / f"{sample}_1.pmta_marker.blastx.tsv"
        b2 = Path(args.blastx_dir) / f"{sample}_2.pmta_marker.blastx.tsv"

        sol_rel = sol_trim = sol_cov = sol_reads = ""
        top_genome = top_genome_rel = top_genome_reads = ""
        n_genomes = 0
        if coverm_path.exists() and coverm_path.stat().st_size > 0:
            coverm = read_table(coverm_path)
            if coverm:
                fields = coverm[0].keys()
                rel_col = find_col(fields, "Relative Abundance")
                trim_col = find_col(fields, "Trimmed Mean")
                cov_col = find_col(fields, "Covered Fraction")
                count_col = find_col(fields, "Read Count")
                vals = []
                for row in coverm:
                    if row.get("Genome") == "unmapped":
                        continue
                    n_genomes += 1
                    rel = float(row.get(rel_col, 0) or 0)
                    trim = float(row.get(trim_col, 0) or 0)
                    cov = float(row.get(cov_col, 0) or 0)
                    count = int(float(row.get(count_col, 0) or 0))
                    vals.append((row.get("Genome"), rel, trim, cov, count))
                if vals:
                    sol_rel = sum(v[1] for v in vals)
                    sol_trim = sum(v[2] for v in vals)
                    sol_cov = sum(v[3] for v in vals) / len(vals)
                    sol_reads = sum(v[4] for v in vals)
                    top = max(vals, key=lambda x: x[1])
                    top_genome, top_genome_rel, top_genome_reads = top[0], top[1], top[4]

        total_reads = ""
        r1_count = count_fastq_reads(r1)
        r2_count = count_fastq_reads(r2)
        if r1_count != "" and r2_count != "":
            total_reads = r1_count + r2_count

        counts = {name: 0 for name in marker_len}
        sample_best = {}
        if total_reads != "":
            for mate, hit_file in ((1, b1), (2, b2)):
                for read_id, hit in best_filtered_hits(hit_file, args.min_identity, args.min_alignment_aa).items():
                    sample_best[f"{mate}:{read_id}"] = hit
            for hit in sample_best.values():
                counts[hit["sseqid"]] = counts.get(hit["sseqid"], 0) + 1
        pmta_hits = sum(counts.values()) if total_reads != "" else ""
        pmta_rpm = pmta_hits / (total_reads / 1_000_000) if total_reads else ""
        pmta_rpkm = pmta_hits / (marker_aa_kb * (total_reads / 1_000_000)) if marker_aa_kb and total_reads else ""
        n_markers_hit = sum(1 for v in counts.values() if v > 0) if total_reads != "" else ""
        top_marker = top_marker_hits = ""
        if counts:
            top_marker, top_marker_hits = max(counts.items(), key=lambda item: item[1])

        detected = ""
        high_conf = ""
        if sol_reads != "" and sol_cov != "":
            detected = "yes" if float(sol_reads) >= 100 and float(sol_cov) >= 0.01 else "no"
            high_conf = "yes" if float(sol_reads) >= 1000 and float(sol_cov) >= 0.05 else "no"

        meta = metadata.get(sample, {})
        rows_out.append(
            {
                "sample_id": sample,
                "record_type": meta.get("record_type", meta.get("source_type", "")),
                "bioproject": meta.get("bioproject", ""),
                "biosample": meta.get("biosample", ""),
                "ecosystem": meta.get("ecosystem", ""),
                "location": meta.get("location", ""),
                "latitude": meta.get("latitude", ""),
                "longitude": meta.get("longitude", ""),
                "pH": meta.get("pH", ""),
                "pH_source": meta.get("pH_source", ""),
                "soil_source_evidence": meta.get("soil_source_evidence", ""),
                "solibacter_relative_abundance_percent": sol_rel,
                "solibacter_mapped_reads": sol_reads,
                "solibacter_trimmed_mean": sol_trim,
                "solibacter_mean_covered_fraction": sol_cov,
                "solibacter_detected": detected,
                "solibacter_high_confidence_detected": high_conf,
                "n_solibacter_genomes": n_genomes,
                "top_solibacter_genome": top_genome,
                "top_solibacter_relative_abundance_percent": top_genome_rel,
                "top_solibacter_read_count": top_genome_reads,
                "pmta_total_reads_searched": total_reads,
                "pmta_marker_hit_reads": pmta_hits,
                "pmta_marker_reads_per_million": pmta_rpm,
                "pmta_marker_rpkm_like": pmta_rpkm,
                "pmta_markers_in_database": len(marker_len),
                "pmta_markers_with_hits": n_markers_hit,
                "pmta_top_marker": top_marker,
                "pmta_top_marker_hits": top_marker_hits,
            }
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"Wrote {len(rows_out)} rows to {out}")


if __name__ == "__main__":
    main()
