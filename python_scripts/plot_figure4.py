import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import random
import yaml

from analysis import run_wilcoxon_diff_abundance, run_maaslin_diff_abundance
from copy import deepcopy
from matplotlib import colors as mcolors
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from utils import trim_taxonomy

with open('config.yml') as f:
    config = yaml.safe_load(f)

def make_volcano_plot(inoculum_signature,
                      drought_signature,
                      num_samples,
                      save_as):
    data = []
    for taxon, vals in inoculum_signature.items():
        pval = vals["adj.P.Val"]
        logfc = vals["logFC"]
        if taxon in drought_signature:
            color = config["colors"]["Treatment"]["Drought"] if drought_signature[taxon]["logFC"] > 0 else config["colors"]["Treatment"]["Control"]
            alpha = 0.9
        else:
            color = "gray"
            alpha = 0.9
        data.append((taxon, logfc, -np.log10(pval), color, alpha))

    df = pd.DataFrame(data, columns=["Taxon", "logFC", "-log10P", "Color", "Alpha"])
    
    df.loc[~((df["logFC"].abs() > 1) & (df["-log10P"] > -np.log10(0.1))), "Alpha"] = 0.6

    plt.figure(figsize=(8,6))
    plt.scatter(df["logFC"], df["-log10P"], c=df["Color"], alpha=df["Alpha"])

    plt.axhline(-np.log10(0.1), color="gray", linestyle="--")
    plt.axvline(-1, color="gray", linestyle="--")
    plt.axvline(1, color="gray", linestyle="--")

    plt.title(f"Samples: {num_samples}")

    plt.xlabel("log2 Fold Change (Inoculum)")
    plt.ylabel("-log10 Adjusted P-Value")
    plt.xlim((-7, 7))
    plt.ylim((0, 15.5))

    # label taxa beyond thresholds
    sig_df = df[(df["logFC"].abs() > 1) & (df["-log10P"] > -np.log10(0.1))]
    for _, row in sig_df.iterrows():
        if row["Taxon"] in drought_signature:
            plt.text(row["logFC"], row["-log10P"],
                     trim_taxonomy(row["Taxon"]),
                     fontsize=8, ha="right", va="bottom")

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()

def permutation_test(drought_signature, inoculum_signature, n_perm=1000, seed=42):
    random.seed(seed)
    np.random.seed(seed)
    
    overlap_taxa = list(set(drought_signature).intersection(inoculum_signature))
    d_signs = np.array([np.sign(drought_signature[t]["logFC"]) for t in overlap_taxa])
    i_signs = np.array([np.sign(inoculum_signature[t]["logFC"]) for t in overlap_taxa])
    
    # observed reversals
    obs = np.sum(d_signs != i_signs)
    
    # null distribution
    null = []
    for _ in range(n_perm):
        permuted = np.random.permutation(d_signs)
        null.append(np.sum(permuted != i_signs))
    
    null = np.array(null)
    
    # one-sided p-value: how often null ≥ observed?
    pval = (np.sum(null >= obs) + 1) / (n_perm + 1)
    
    return obs, len(overlap_taxa), pval

def main():
    success_studies = ["moore2023microbial", "zhang2022cross"]
    fail_studies = ["munoz-ucros2021drought", "swift2024drought"]

    all_studies = success_studies + fail_studies

    merged_signature = {}
    # Get compartment and associated signature
    for compartment in ["bulk_soil", "endosphere", "rhizosphere"]:
        maaslin_signature_file = f"../data/maaslin_signature_{compartment}.pkl"
        with open(maaslin_signature_file, 'rb') as handle:
            compartment_signature = pkl.load(handle)
            compartment_signature = {taxon: compartment_signature[taxon]
                                     for taxon in compartment_signature
                                     if compartment_signature[taxon]["adj.P.Val"] <= 0.05}

        for taxon, vals in compartment_signature.items():
            logFC = vals["logFC"]
            adjP = vals["adj.P.Val"]

            if taxon not in merged_signature:
                merged_signature[taxon] = {"logFC": logFC, "adj.P.Val": adjP, "signs": {1 if logFC > 0 else -1}}
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

    print(len(merged_signature))

    inverses = dict()
    for study in all_studies:
        print(study)
        # Load dataset
        pkl_file = f"../data/{study}_l6_filtered.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        # Filter for plant-associated samples
        ds.filter_rows(lambda df: df["IsPlantAssociated"] == "Yes")

        # Filter out sterile controls
        ds.filter_rows(lambda df: df["Inoculum"] != "Sterile")

        # Filter low-abundance features
        ds.filter_features()

        # Run wilcoxon
        inoculum_signature = run_maaslin_diff_abundance(ds,
                                                        fixed_effects="Inoculum",
                                                        random_effects="Treatment",
                                                        base="Control",
                                                        condition="DroughtLegacy")

        inoculum_signature = {taxon.replace(".", ";"): inoculum_signature[taxon]
                              for taxon in inoculum_signature}


        make_volcano_plot(inoculum_signature,
                          merged_signature,
                          num_samples=len(ds.taxonomy_counts_df),
                          save_as=f"../data/plots/volcano_{study}")

        inverses[study] = permutation_test(merged_signature, inoculum_signature)

    print(inverses)
                

if __name__ == "__main__":
    main()
