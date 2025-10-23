import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import re
import yaml

from analysis import run_maaslin_diff_abundance, run_wilcoxon_diff_abundance, run_mmuphin_diff_abundance
from copy import deepcopy
from matplotlib import colors as mcolors
from utils import trim_taxonomy

with open('config.yml') as f:
    config = yaml.safe_load(f)

def plot_signature(signature,
                   save_as):
    data = [(taxon, vals["logFC"]) for taxon, vals in signature.items()]
    df = pd.DataFrame(data, columns=["Taxon", "logFC"])
    df = pd.concat([df[df["logFC"] > 0].sort_values("logFC", ascending=False),
                    df[df["logFC"] < 0].sort_values("logFC", ascending=False)])
    
    colors = [config["colors"]["Treatment"]["Drought"] if v > 0
              else config["colors"]["Treatment"]["Control"] for v in df["logFC"]]
    
    plt.figure(figsize=(0.4*len(df), 6))
    plt.bar(df["Taxon"], df["logFC"], color=colors)

    print(signature)
    trimmed_labels = [trim_taxonomy(t) for t in df["Taxon"]]
    plt.xticks(range(len(df)), trimmed_labels, rotation=90)
    plt.ylim((-2.1, 1.5))

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()

    return df["Taxon"]

def normalize_name(name):
    return re.sub(r'[^A-Za-z0-9]', '', name).lower()

def plot_per_study_validation(study_signatures, main_signature, save_as):
    studies = list(study_signatures.keys())
    taxa = main_signature

    fig, ax = plt.subplots(figsize=(max(0.4*len(taxa), 8), 0.4*len(studies)))

    for i, study in enumerate(studies):
        sig = study_signatures[study]
        # Build a normalized lookup dictionary
        normalized_sig = {normalize_name(k): v for k, v in sig.items()}

        for j, taxon in enumerate(taxa):
            normalized_taxon = normalize_name(taxon)
            if normalized_taxon not in normalized_sig:
                facecolor = "white"
            else:
                entry = normalized_sig[normalized_taxon]
                logfc = entry["logFC"]
                pval = entry["adj.P.Val"]

                if logfc > 0:
                    base_color = config["colors"]["Treatment"]["Drought"]
                elif logfc < 0:
                    base_color = config["colors"]["Treatment"]["Control"]
                else:
                    base_color = "white"

                if pval < 0.05:
                    facecolor = base_color
                else:
                    facecolor = mcolors.to_rgba(base_color, alpha=0.5)

            ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=facecolor, edgecolor="none"))

    ax.set_xticks(np.arange(len(taxa)) + 0.5)
    ax.set_yticks(np.arange(len(studies)) + 0.5)
    ax.set_xticklabels([trim_taxonomy(t) for t in taxa], rotation=90)
    ax.set_yticklabels(studies)

    ax.set_xlim(0, len(taxa))
    ax.set_ylim(0, len(studies))
    ax.invert_yaxis()
    ax.set_aspect("equal")

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def main():
    for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus_for_signature_batch_corrected.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        maaslin_signature_file = f"../raw_data/maaslin_signature_{compartment.lower().replace(' ', '_')}.pkl"
        maaslin_signature = run_mmuphin_diff_abundance(ds)
        with open(maaslin_signature_file , 'wb') as handle:
            pkl.dump(maaslin_signature, handle)

        # Save signature as tsv
        df = pd.DataFrame.from_dict(maaslin_signature,
                                    orient='index')
        df.index.name = "Taxon"
        tsv_file = f"../raw_data/signature_{compartment.lower().replace(' ', '_')}.tsv"
        df.to_csv(tsv_file, sep='\t')

        filtered_maaslin_signature = {taxon: maaslin_signature[taxon]
                                      for taxon in maaslin_signature
                                      if maaslin_signature[taxon]["adj.P.Val"] <= 0.05}

        print(f"Maaslin {len(filtered_maaslin_signature)}")
    
        signature_taxa = plot_signature(filtered_maaslin_signature,
                                        save_as=f"../plots/maaslin_signature_{compartment.lower().replace(' ', '_')}")

        # Per-study validation      
        study_signatures = dict()
        for study in config["all_studies"]:
            pkl_file = f"../datasets/{study}/genus.pkl"
            with open(pkl_file, 'rb') as handle:
                study_ds = pkl.load(handle)

            study_ds.filter_samples_based_on_condition(lambda df: df["RootCompartment"] == compartment)
            study_ds.remove_zero_features()

            if not len(study_ds.taxonomy_counts_df):
                continue
                
            study_signature = run_maaslin_diff_abundance(study_ds)
            # For some reason - and ; are replaced by .
            study_signatures[study] = {key.replace(".", ";"): study_signature[key]
                                       for key in study_signature}

        print(study_signatures)
        plot_per_study_validation(study_signatures,
                                  signature_taxa,
                                  save_as=f"../plots/signature_validation_{compartment.lower().replace(' ', '_')}")

        # For the supplement: PERMANOVA results before and after batch effect correction
        # process_permanova(6, compartment.lower().replace(' ', '_'), config, True)
        # process_permanova(6, compartment.lower().replace(' ', '_'), config, False)

if __name__ == "__main__":
    main()
