from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGETS = ROOT / "config" / "targets.fna"
OUT_DIR = ROOT / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                records[name] = "".join(chunks).upper()
            name = line[1:].split()[0]
            chunks = []
        else:
            chunks.append(line)
    if name:
        records[name] = "".join(chunks).upper()
    return records


def translate(seq):
    protein = []
    for i in range(0, len(seq) - 2, 3):
        protein.append(CODON_TABLE.get(seq[i : i + 3], "X"))
    return "".join(protein)


def wrap(seq, width=60):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


records = read_fasta(TARGETS)
pmtA_nt = records["Ellin6076_pmtA_seed"]
pmtA_aa = translate(pmtA_nt)

(OUT_DIR / "Ellin6076_pmtA_seed.fna").write_text(
    f">Ellin6076_pmtA_seed\n{wrap(pmtA_nt)}\n", encoding="utf-8"
)
(OUT_DIR / "Ellin6076_pmtA_seed.faa").write_text(
    f">Ellin6076_pmtA_seed_translated\n{wrap(pmtA_aa)}\n", encoding="utf-8"
)
(OUT_DIR / "Ellin6076_16S_partial.fna").write_text(
    f">Ellin6076_16S_partial\n{wrap(records['Ellin6076_16S_partial'])}\n",
    encoding="utf-8",
)

summary = [
    "target\tlength_nt\tlength_aa\tstarts_with_atg\tterminal_stop\tinternal_stop_count",
    "\t".join(
        [
            "Ellin6076_pmtA_seed",
            str(len(pmtA_nt)),
            str(len(pmtA_aa.rstrip("*"))),
            str(pmtA_nt.startswith("ATG")),
            str(pmtA_aa.endswith("*")),
            str(pmtA_aa[:-1].count("*")),
        ]
    ),
    "\t".join(
        [
            "Ellin6076_16S_partial",
            str(len(records["Ellin6076_16S_partial"])),
            "NA",
            "NA",
            "NA",
            "NA",
        ]
    ),
]
(OUT_DIR / "target_sequence_qc.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")
print("\n".join(summary))
