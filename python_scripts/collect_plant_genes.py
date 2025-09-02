import os
import csv
from pathlib import Path

monocots = {
    "Oryza sativa", "Echinochloa esculenta", "Setaria italica", "Sorghum bicolor", "Panicum miliaceum",
    "Cenchrus americanus", "Lolium arundinaceum", "Phalaris arundinacea", "Sporobolus cryptandrus",
    "Bromus inermis", "Phleum pratense", "Sporobolus neglectus", "Digitaria ischaemum", "Oryza glaberrima",
    "Triticum aestivum", "Bothriochloa bladhii", "Sorghastrum nutans", "Brachypodium distachyon", "Avena sativa",
    "Triticum monococcum", "Secale cereale", "Triticum turgidum", "Hordeum vulgare", "Zea mays", "Eragrostis tef",
    "Miscanthus sinensis", "Panicum virgatum"
}

def load_annotations(tsv_path):
    annotations = {}
    with open(tsv_path, newline='') as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            plant = row["Plant"].replace(" ", "_")
            chloro = row["Chloroplast genome"].strip()
            rbcl_acc = chloro if chloro else row["rbcL"].strip()
            matk_acc = chloro if chloro else row["matK"].strip()
            its_acc = row["ITS1_5.8S_ITS2"].strip()
            annotations[plant] = {
                "rbcl": rbcl_acc if rbcl_acc else "NA",
                "matk": matk_acc if matk_acc else "NA",
                "its": its_acc if its_acc else "NA",
            }
    return annotations

def collect_sequences_split(locus, output_file_prefix, annotations):
    monocot_file = f"{output_file_prefix}_monocots.fasta"
    dicot_file = f"{output_file_prefix}_dicots.fasta"

    with open(monocot_file, "w") as mono_f, open(dicot_file, "w") as dicot_f:
        for plant_dir in Path("../data/plant_genes").iterdir():
            plant = plant_dir.name
            seq_file = plant_dir / f"{locus}.fna"
            if seq_file.exists() and plant in annotations:
                accession = annotations[plant][locus] if annotations[plant][locus] != "NA" else "unknown"
                with open(seq_file) as f:
                    lines = f.readlines()
                    seq = "".join(lines[1:]).replace("\n", "")
                    entry = f">{plant}|{locus}|{accession}\n{seq}\n"
                    if plant.replace("_", " ") in monocots:
                        mono_f.write(entry)
                    else:
                        dicot_f.write(entry)

annotations = load_annotations("../data/plant_gene_accessions.tsv")
collect_sequences_split("rbcl", "../data/plant_genes/all_rbcl", annotations)
collect_sequences_split("matk", "../data/plant_genes/all_matk", annotations)
collect_sequences_split("its", "../data/plant_genes/all_its", annotations)
