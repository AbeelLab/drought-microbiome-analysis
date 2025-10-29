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

# --- Utility functions ---
def normalize_name(name):
    return re.sub(r'[^A-Za-z0-9]', '', name).lower()


# --- Signature subplot function ---
def plot_all_signatures(signatures_dict, save_as):
    compartments = list(signatures_dict.keys())
    n_compartments = len(compartments)

    fig, axes = plt.subplots(1, n_compartments, figsize=(6 * n_compartments, 6), sharey=True)

    if n_compartments == 1:
        axes = [axes]

    for ax, compartment in zip(axes, compartments):
        signature = signatures_dict[compartment]
        if not signature:
            continue

        data = [(taxon, vals["logFC"]) for taxon, vals in signature.items()]
        df = pd.DataFrame(data, columns=["Taxon", "logFC"])
        df = pd.concat([
            df[df["logFC"] > 0].sort_values("logFC", ascending=False),
            df[df["logFC"] < 0].sort_values("logFC", ascending=False)
        ])

        colors = [config["colors"]["Treatment"]["Drought"] if v > 0
                  else config["colors"]["Treatment"]["Control"] for v in df["logFC"]]

        ax.bar(df["Taxon"], df["logFC"], color=colors)
        trimmed_labels = [trim_taxonomy(t) for t in df["Taxon"]]
        ax.set_xticks(range(len(df)))
        ax.set_xticklabels(trimmed_labels, rotation=90)
        ax.set_ylim((-2.1, 1.5))
        ax.set_title(compartment)

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


# --- Validation subplot function ---
def plot_all_study_validations(all_study_signatures, main_signatures, save_as):
    compartments = list(main_signatures.keys())
    all_studies = config["all_studies"]

    fig, axes = plt.subplots(1, len(compartments), figsize=(8 * len(compartments), 0.6 * len(all_studies)))

    if len(compartments) == 1:
        axes = [axes]

    for ax, compartment in zip(axes, compartments):
        study_signatures = all_study_signatures.get(compartment, {})
        taxa = list(main_signatures[compartment].keys())

        for i, study in enumerate(all_studies):
            if study not in study_signatures:
                # Gray squares (no data for this compartment)
                for j in range(len(taxa)):
                    ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="gray", edgecolor="none"))
                continue

            sig = study_signatures[study]
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

                    facecolor = base_color if pval < 0.05 else mcolors.to_rgba(base_color, alpha=0.5)

                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=facecolor, edgecolor="none"))

        ax.set_xticks(np.arange(len(taxa)) + 0.5)
        ax.set_xticklabels([trim_taxonomy(t) for t in taxa], rotation=90)
        ax.set_yticks(np.arange(len(all_studies)) + 0.5)
        ax.set_yticklabels(all_studies)
        ax.set_xlim(0, len(taxa))
        ax.set_ylim(0, len(all_studies))
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set_title(compartment)

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


# --- Main pipeline ---
def main():
    all_signatures = {}
    all_study_signatures = {}

    for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus_for_signature_batch_corrected.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        maaslin_signature_file = f"../raw_data/maaslin_signature_{compartment.lower().replace(' ', '_')}.pkl"
        maaslin_signature = run_mmuphin_diff_abundance(ds)
        with open(maaslin_signature_file, 'wb') as handle:
            pkl.dump(maaslin_signature, handle)

        # Save signature as TSV
        df = pd.DataFrame.from_dict(maaslin_signature, orient='index')
        df.index.name = "Taxon"
        tsv_file = f"../raw_data/signature_{compartment.lower().replace(' ', '_')}.tsv"
        df.to_csv(tsv_file, sep='\t')

        filtered_maaslin_signature = {
            taxon: maaslin_signature[taxon]
            for taxon in maaslin_signature
            if maaslin_signature[taxon]["adj.P.Val"] <= 0.05
        }

        print(f"{compartment} significant taxa: {len(filtered_maaslin_signature)}")

        all_signatures[compartment] = filtered_maaslin_signature

        # Per-study validation
        study_signatures = {}
        for study in config["all_studies"]:
            pkl_file = f"../datasets/{study}/genus.pkl"
            if not os.path.exists(pkl_file):
                continue

            with open(pkl_file, 'rb') as handle:
                study_ds = pkl.load(handle)

            study_ds.filter_samples_based_on_condition(lambda df: df["RootCompartment"] == compartment)
            study_ds.remove_zero_features()

            if not len(study_ds.taxonomy_counts_df):
                continue

            study_signature = run_maaslin_diff_abundance(study_ds)
            study_signatures[study] = {
                key.replace(".", ";"): study_signature[key]
                for key in study_signature
            }

        all_study_signatures[compartment] = study_signatures

    # --- Combined plots ---
    plot_all_signatures(all_signatures, save_as="../plots/maaslin_signatures_all")
    plot_all_study_validations(all_study_signatures, all_signatures, save_as="../plots/signature_validation_all")


if __name__ == "__main__":
    main()

