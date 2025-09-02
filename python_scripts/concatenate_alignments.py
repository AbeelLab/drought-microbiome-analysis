from Bio import SeqIO
from pathlib import Path

def parse_alignment(path):
    path = Path(path)
    if not path.exists():
        return {}
    data = {}
    for record in SeqIO.parse(str(path), "fasta"):
        plant = record.id.split("|")[0]
        data[plant] = str(record.seq)
    return data

def build_supermatrix(out_prefix, rbcl_file, matk_file, its_file):
    rbcl = parse_alignment(rbcl_file)
    matk = parse_alignment(matk_file)
    its = parse_alignment(its_file)

    rbcl_len = len(next(iter(rbcl.values()))) if rbcl else 0
    matk_len = len(next(iter(matk.values()))) if matk else 0
    its_len = len(next(iter(its.values()))) if its else 0

    all_plants = sorted(set(rbcl) | set(matk) | set(its))

    out_path = Path(out_prefix).with_suffix(".fasta")
    with open(out_path, "w") as out_f:
        for plant in all_plants:
            seq_rbcl = rbcl.get(plant, "-" * rbcl_len)
            seq_matk = matk.get(plant, "-" * matk_len)
            seq_its = its.get(plant, "-" * its_len)
            full_seq = seq_rbcl + seq_matk + seq_its
            out_f.write(f">{plant}\n{full_seq}\n")

base = Path("../data/plant_genes")

build_supermatrix(
    base / "all_aligned_genes",
    base / "aligned_rbcl.fasta",
    base / "aligned_matk.fasta",
    base / "aligned_its.fasta",
)

build_supermatrix(
    base / "all_aligned_genes_monocots",
    base / "aligned_rbcl_monocots.fasta",
    base / "aligned_matk_monocots.fasta",
    base / "aligned_its_monocots.fasta",
)

build_supermatrix(
    base / "all_aligned_genes_dicots",
    base / "aligned_rbcl_dicots.fasta",
    base / "aligned_matk_dicots.fasta",
    base / "aligned_its_dicots.fasta",
)
