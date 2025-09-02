import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import yaml

from analysis import run_wilcoxon_diff_abundance, run_mmuphin_diff_abundance, run_limma_diff_abundance
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

    trimmed_labels = [trim_taxonomy(t) for t in df["Taxon"]]
    plt.xticks(range(len(df)), trimmed_labels, rotation=90)

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()

    return df["Taxon"]

def plot_per_study_validation(study_signatures, main_signature, save_as):
    studies = list(study_signatures.keys())
    taxa = main_signature

    fig, ax = plt.subplots(figsize=(max(0.4*len(taxa), 8), 0.4*len(studies)))

    for i, study in enumerate(studies):
        sig = study_signatures[study]
        for j, taxon in enumerate(taxa):
            if taxon not in sig:
                facecolor = "white"
                print("***", taxon)
            else:
                logfc = sig[taxon]["logFC"]
                pval = sig[taxon]["adj.P.Val"]
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
        pkl_file = f"../data/merged_l6_{compartment.lower().replace(' ', '_')}_intersection_batch_corrected.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        maaslin_signature_file = f"../data/maaslin_signature_{compartment.lower().replace(' ', '_')}.pkl"
        maaslin_signature = run_mmuphin_diff_abundance(ds)
        with open(maaslin_signature_file , 'wb') as handle:
            pkl.dump(maaslin_signature, handle)

        if compartment != "Bulk soil":
            formula = "Treatment + StudyID + HostSpecific"
        else:
            formula = "Treatment + StudyID"
        limma_signature_file = f"../data/limma_signature_{compartment.lower().replace(' ', '_')}.pkl"
        limma_signature = run_limma_diff_abundance(ds,
                                                   formula=formula)
        with open(limma_signature_file , 'wb') as handle:
            pkl.dump(limma_signature, handle)

        filtered_limma_signature = {taxon: limma_signature[taxon]
                                    for taxon in limma_signature
                                    if limma_signature[taxon]["adj.P.Val"] <= 0.05}
        filtered_maaslin_signature = {taxon: maaslin_signature[taxon]
                                      for taxon in maaslin_signature
                                      if maaslin_signature[taxon]["adj.P.Val"] <= 0.05}

        print(f"Maaslin {len(filtered_maaslin_signature)}")
        print(f"Limma {len(filtered_limma_signature)}")
        print(f"{compartment} overlap: {len(set(filtered_maaslin_signature.keys()).intersection(set(filtered_limma_signature.keys())))}")
    
        signature_taxa = plot_signature(filtered_maaslin_signature,
                                        save_as=f"../data/plots/maaslin_signature_{compartment.lower().replace(' ', '_')}")

        # Per-study validation
        # Load non-batch-corrected dataset
        pkl_file = f"../data/merged_l6_unfiltered_{compartment.lower().replace(' ', '_')}_intersection.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)
        study_signatures = dict()
        for study in ds.metadata_df["StudyName"].unique().tolist():
            study_ds = deepcopy(ds)
            study_ds.filter_rows(lambda df: df["StudyName"] == study)
            study_ds.remove_zero_features()

            study_signature = run_wilcoxon_diff_abundance(study_ds)
            study_signatures[study] = study_signature
        plot_per_study_validation(study_signatures,
                                  signature_taxa,
                                  save_as=f"../data/plots/signature_validation_{compartment.lower().replace(' ', '_')}")

if __name__ == "__main__":
    main()
