#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


SAM_RE = re.compile(r"[ED].G.G.G")
HYDROPHOBIC = set("AILMFWVY")
BASIC = set("KRH")


def read_fasta(path):
    records = {}
    name = None
    seqs = []
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seqs).replace("*", "")
                name = line[1:].split()[0]
                seqs = []
            else:
                seqs.append(line)
    if name is not None:
        records[name] = "".join(seqs).replace("*", "")
    return records


def wrap(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def nterm_membrane_proxy(seq, window=35):
    """Heuristic proxy for an N-terminal amphipathic/membrane-binding helix.

    This is not a structural prediction. It flags proteins whose N terminus is
    enriched for hydrophobic residues and contains basic residues, consistent
    with the conserved PmtA N-terminal membrane-binding feature.
    """
    nterm = seq[:window]
    if len(nterm) < 20:
        return False, 0.0, 0.0
    hydrophobic_fraction = sum(aa in HYDROPHOBIC for aa in nterm) / len(nterm)
    basic_fraction = sum(aa in BASIC for aa in nterm) / len(nterm)
    passes = hydrophobic_fraction >= 0.35 and basic_fraction >= 0.03
    return passes, hydrophobic_fraction, basic_fraction


def motif_positions(seq):
    return [m.start() + 1 for m in SAM_RE.finditer(seq)]


def classify_hit(pident, qcov, scov, evalue, has_sam, has_nterm):
    if pident >= 30 and qcov >= 0.60 and evalue <= 1e-10 and has_sam and has_nterm:
        return "high_confidence_functional_PmtA_like"
    if pident >= 25 and qcov >= 0.45 and evalue <= 1e-5 and has_sam:
        return "PmtA_like_requires_manual_review"
    return "reject"


parser = argparse.ArgumentParser(
    description="Filter gene101949-seeded PmtA protein hits by conserved protein features."
)
parser.add_argument("--blast", required=True, help="DIAMOND BLASTP outfmt 6 table.")
parser.add_argument("--proteins", required=True, help="Predicted protein FASTA used as DIAMOND database.")
parser.add_argument("--functional-seed", required=True, help="Functional gene101949 PmtA seed FASTA.")
parser.add_argument("--paralog-seed", required=True, help="gene102575 PmtA-like paralog FASTA for reference output.")
parser.add_argument("--out-tsv", required=True)
parser.add_argument("--out-faa", required=True)
args = parser.parse_args()

proteins = read_fasta(args.proteins)
functional_seed = read_fasta(args.functional_seed)
paralog_seed = read_fasta(args.paralog_seed)

fieldnames = [
    "qseqid",
    "sseqid",
    "pident",
    "alignment_length",
    "qcov",
    "scov",
    "evalue",
    "bitscore",
    "slen",
    "has_sam_binding_motif_EDXGXGXG",
    "sam_motif_positions",
    "has_nterm_membrane_proxy",
    "nterm_hydrophobic_fraction",
    "nterm_basic_fraction",
    "classification",
]

best_by_subject = {}
with Path(args.blast).open(encoding="utf-8") as fh:
    for line in fh:
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        (
            qseqid,
            sseqid,
            pident,
            length,
            _mismatch,
            _gapopen,
            _qstart,
            _qend,
            _sstart,
            _send,
            evalue,
            bitscore,
            qlen,
            slen,
        ) = parts
        pident_f = float(pident)
        length_i = int(length)
        qlen_i = int(qlen)
        slen_i = int(slen)
        evalue_f = float(evalue)
        bitscore_f = float(bitscore)
        qcov = length_i / qlen_i if qlen_i else 0.0
        scov = length_i / slen_i if slen_i else 0.0
        current = best_by_subject.get(sseqid)
        if current is None or bitscore_f > current["bitscore"]:
            best_by_subject[sseqid] = {
                "qseqid": qseqid,
                "sseqid": sseqid,
                "pident": pident_f,
                "alignment_length": length_i,
                "qcov": qcov,
                "scov": scov,
                "evalue": evalue_f,
                "bitscore": bitscore_f,
                "slen": slen_i,
            }

rows = []
for sseqid, hit in sorted(best_by_subject.items()):
    seq = proteins.get(sseqid, "")
    sam_positions = motif_positions(seq)
    has_sam = bool(sam_positions)
    has_nterm, hydrophobic_fraction, basic_fraction = nterm_membrane_proxy(seq)
    classification = classify_hit(
        hit["pident"], hit["qcov"], hit["scov"], hit["evalue"], has_sam, has_nterm
    )
    row = dict(hit)
    row.update(
        {
            "has_sam_binding_motif_EDXGXGXG": has_sam,
            "sam_motif_positions": ";".join(map(str, sam_positions)),
            "has_nterm_membrane_proxy": has_nterm,
            "nterm_hydrophobic_fraction": hydrophobic_fraction,
            "nterm_basic_fraction": basic_fraction,
            "classification": classification,
        }
    )
    rows.append(row)

out_tsv = Path(args.out_tsv)
out_tsv.parent.mkdir(parents=True, exist_ok=True)
with out_tsv.open("w", encoding="utf-8", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)

out_faa = Path(args.out_faa)
out_faa.parent.mkdir(parents=True, exist_ok=True)
with out_faa.open("w", encoding="utf-8") as fh:
    for name, seq in functional_seed.items():
        fh.write(f">{name}|seed=functional_gene101949\n{wrap(seq)}\n")
    for row in rows:
        if row["classification"] == "high_confidence_functional_PmtA_like":
            seq = proteins.get(row["sseqid"], "")
            if seq:
                fh.write(f">{row['sseqid']}|classification={row['classification']}\n{wrap(seq)}\n")

print(f"Wrote {len(rows)} PmtA-like candidate rows to {out_tsv}")
print(f"Wrote functional marker FASTA to {out_faa}")
print("Note: gene102575 is treated as a PmtA-like paralog/control, not as the functional seed.")
