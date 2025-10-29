import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import yaml

from analysis import run_wilcoxon_diff_abundance, run_maaslin_diff_abundance
from copy import deepcopy
from matplotlib import colors as mcolors
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from utils import trim_taxonomy

# Load configuration
with open('config.yml') as f:
    config = yaml.safe_load(f)


def make_volcano_subplots(inoculum_signatures,
                          drought_signature,
                          inoculum_taxa,
                          processed_studies,
                          save_as):
    """Generate a combined multi-panel volcano plot for all studies."""
    studies = list(inoculum_signatures.keys())
    n_studies = len(studies)
    ncols = 2
    nrows = int(np.ceil(n_studies / ncols))

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12, 12))
    axes = axes.flatten()

    for i, study in enumerate(studies):
        ax = axes[i]
        inoculum_signature = inoculum_signatures[study]
        num_samples = len(processed_studies[study].taxonomy_counts_df)

        data = []
        for taxon, vals in inoculum_signature.items():
            pval = vals["adj.P.Val"]
            logfc = vals["logFC"]
            if taxon in drought_signature:
                color = (
                    config["colors"]["Treatment"]["Drought"]
                    if drought_signature[taxon]["logFC"] > 0
                    else config["colors"]["Treatment"]["Control"]
                )
                alpha = 0.9
            elif taxon in inoculum_taxa:
                color = "black"
                alpha = 0.9
            else:
                color = "gray"
                alpha = 0.6
            data.append((taxon, logfc, -np.log10(pval), color, alpha))

        df = pd.DataFrame(data, columns=["Taxon", "logFC", "-log10P", "Color", "Alpha"])
        df.loc[~((df["logFC"].abs() > 1.75) & (df["-log10P"] > -np.log10(0.05))), "Alpha"] = 0.5

        ax.scatter(df["logFC"], df["-log10P"], c=df["Color"], alpha=df["Alpha"], s=20)

        ax.axhline(-np.log10(0.05), color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(-1.75, color="gray", linestyle="--", linewidth=0.8)
        ax.axvline(1.75, color="gray", linestyle="--", linewidth=0.8)

        ax.set_title(f"{study}\nSamples: {num_samples}", fontsize=12)
        ax.set_xlabel("log2 Fold Change (Inoculum)", fontsize=11)
        ax.set_ylabel("-log10 Adjusted P-Value", fontsize=11)
        ax.tick_params(axis="both", labelsize=10)

        ax.set_xlim((-7, 7))
        ax.set_ylim((0, 15.5))

        # Label significant taxa
        sig_df = df[(df["logFC"].abs() > 1.75) & (df["-log10P"] > -np.log10(0.05))]
        for _, row in sig_df.iterrows():
            if row["Taxon"] in drought_signature or row["Taxon"] in inoculum_taxa:
                ax.text(row["logFC"], row["-log10P"],
                        trim_taxonomy(row["Taxon"]),
                        fontsize=7, ha="right", va="bottom")

    # Hide unused subplots
    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.tight_layout()
    fig.savefig(f"{save_as}.svg", format="svg")
    fig.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close(fig)


def permutation_test(drought_signature, inoculum_signature, n_perm=100):
    overlap_taxa = list(set(drought_signature).intersection(inoculum_signature))
    d_signs = np.array([np.sign(drought_signature[t]["logFC"]) for t in overlap_taxa])
    i_signs = np.array([np.sign(inoculum_signature[t]["logFC"]) for t in overlap_taxa])

    matrix = {
        "drought_depleted_inoculum_depleted": np.sum((d_signs == -1) & (i_signs == -1)),
        "drought_depleted_inoculum_enriched": np.sum((d_signs == -1) & (i_signs == 1)),
        "drought_enriched_inoculum_depleted": np.sum((d_signs == 1) & (i_signs == -1)),
        "drought_enriched_inoculum_enriched": np.sum((d_signs == 1) & (i_signs == 1)),
    }

    print(f"Overlap taxa: {len(overlap_taxa)}")
    print("Contingency matrix (counts):")
    print(f"{'':20s}{'Inoculum depleted':>20s}{'Inoculum enriched':>20s}")
    print(f"{'Drought depleted':20s}{matrix['drought_depleted_inoculum_depleted']:>20}{matrix['drought_depleted_inoculum_enriched']:>20}")
    print(f"{'Drought enriched':20s}{matrix['drought_enriched_inoculum_depleted']:>20}{matrix['drought_enriched_inoculum_enriched']:>20}")

    obs = np.sum(d_signs != i_signs)
    null = []
    for i in range(n_perm):
        np.random.seed(i)
        permuted = np.random.permutation(d_signs)
        null.append(np.sum(permuted != i_signs))
    null = np.array(null)
    pval = (np.sum(null >= obs) + 1) / (n_perm + 1)

    print(f"Inversions p-val: {pval}")

    return obs, len(overlap_taxa), pval


def main():
    success_studies = ["moore2023microbial", "zhang2022cross"]
    fail_studies = ["munoz-ucros2022drought", "swift2025drought"]

    all_studies = success_studies + fail_studies

    merged_signature = {}
    # Merge compartment signatures
    for compartment in ["bulk_soil", "endosphere", "rhizosphere"]:
        maaslin_signature_file = f"../raw_data/maaslin_signature_{compartment}.pkl"
        with open(maaslin_signature_file, 'rb') as handle:
            compartment_signature = pkl.load(handle)
            compartment_signature = {
                taxon: compartment_signature[taxon]
                for taxon in compartment_signature
                if compartment_signature[taxon]["adj.P.Val"] <= 0.05
            }

        for taxon, vals in compartment_signature.items():
            logFC = vals["logFC"]
            adjP = vals["adj.P.Val"]

            if taxon not in merged_signature:
                merged_signature[taxon] = {
                    "logFC": logFC,
                    "adj.P.Val": adjP,
                    "signs": {1 if logFC > 0 else -1}
                }
            else:
                current = merged_signature[taxon]
                current["signs"].add(1 if logFC > 0 else -1)
                if abs(logFC) > abs(current["logFC"]):
                    current["logFC"] = logFC
                current["adj.P.Val"] = max(current["adj.P.Val"], adjP)

    for taxon, vals in merged_signature.items():
        if len(vals["signs"]) > 1:
            vals["logFC"] = 0
        del vals["signs"]

    inoculum_signatures = {}
    inoculum_signature_counts = {}
    processed_studies = {}

    for study in all_studies:
        print(f"[INFO] Processing study: {study}")
        pkl_file = f"../datasets/{study}/genus.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        ds.filter_samples_based_on_condition(lambda df: df["Inoculum"] != "Sterile")
        ds.filter_features()
        processed_studies[study] = ds

        inoculum_signature = run_maaslin_diff_abundance(
            ds,
            fixed_effects="Inoculum",
            random_effects="Treatment",
            base="Reference",
            condition="Dry"
        )

        inoculum_signature = {
            taxon.replace(".", ";"): inoculum_signature[taxon]
            for taxon in inoculum_signature
        }

        inoculum_signatures[study] = inoculum_signature

        filtered_inoculum_signature = {
            taxon: inoculum_signature[taxon]
            for taxon in inoculum_signature
            if inoculum_signature[taxon]["adj.P.Val"] <= 0.05 and abs(inoculum_signature[taxon]["logFC"]) >= 1.75
        }

        for taxon in filtered_inoculum_signature:
            inoculum_signature_counts[taxon] = inoculum_signature_counts.get(taxon, 0) + 1

    filtered_inoculum_taxa = {taxon for taxon, count in inoculum_signature_counts.items() if count > 1}

    # ---- Multi-panel plot ----
    make_volcano_subplots(
        inoculum_signatures=inoculum_signatures,
        drought_signature=merged_signature,
        inoculum_taxa=filtered_inoculum_taxa,
        processed_studies=processed_studies,
        save_as="../plots/volcano_all_studies"
    )

    # ---- Permutation test for each study ----
    for study in all_studies:
        permutation_test(merged_signature, inoculum_signatures[study])


if __name__ == "__main__":
    main()
