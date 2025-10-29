import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle as pkl
import seaborn as sns
import yaml

from copy import deepcopy
from matplotlib.patches import Patch
from matplotlib_venn import venn3
from sklearn.metrics import jaccard_score

with open('config.yml') as f:
    config = yaml.safe_load(f)

def plot_jaccard_heatmap(ds, save_as, figsize=(8, 6), cmap="Greys"):
    metadata = ds.metadata_df.copy()
    study_names = metadata["StudyName"]

    presence_absence = (ds.taxonomy_counts_df > 0).astype(int)

    study_pa = presence_absence.groupby(study_names).max()

    studies = study_pa.index
    n = len(studies)
    jaccard_matrix = np.zeros((n, n))

    for i, s1 in enumerate(studies):
        for j, s2 in enumerate(studies):
            if i <= j:
                score = jaccard_score(study_pa.loc[s1], study_pa.loc[s2])
                jaccard_matrix[i, j] = score
                jaccard_matrix[j, i] = score

    jaccard_df = pd.DataFrame(jaccard_matrix, index=studies, columns=studies)

    union_features = (study_pa.sum(axis=0) > 0).sum()
    intersection_features = (study_pa.sum(axis=0) == len(studies)).sum()

    plt.figure(figsize=figsize)
    sns.heatmap(
        jaccard_df,
        annot=True,
        fmt=".2f",
        cmap=cmap,
        vmin=0,
        vmax=1,
        cbar_kws={"label": "Jaccard Similarity"}
    )

    plt.title("Jaccard Similarity Between Studies", pad=20)

    plt.figtext(
        0.5, -0.05,
        f"Genera prevalent in at least one study: {union_features}\n"
        f"Genera prevalent in all studies: {intersection_features}",
        ha="center", va="top", fontsize=10)
    print("Union features", union_features)
    print("Intersection features", intersection_features)

    plt.tight_layout()

    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()

    return jaccard_df, union_features, intersection_features

def plot_rarefaction_per_compartment(save_as, n_permutations=100, figsize=(7, 8)):
    plt.figure(figsize=figsize)

    np.random.seed(42)
    for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        counts = ds.taxonomy_counts_df.copy()
        n_samples = counts.shape[0]
        all_curves = np.zeros((n_permutations, n_samples))

        for perm in range(n_permutations):
            sample_order = np.random.permutation(counts.index)
            seen_features = set()
            curve = []

            for i, sample in enumerate(sample_order):
                nonzero_features = set(counts.loc[sample][counts.loc[sample] > 0].index)
                seen_features.update(nonzero_features)
                curve.append(len(seen_features))

            all_curves[perm, :] = curve

        mean_curve = all_curves.mean(axis=0)
        std_curve = all_curves.std(axis=0)

        x = np.arange(1, n_samples + 1)
        color = config["colors"]["RootCompartment"][compartment]

        plt.plot(x, mean_curve, label=compartment, color=color, linewidth=2.5)
        plt.fill_between(x, mean_curve - std_curve, mean_curve + std_curve,
                         color=color, alpha=0.3)

    plt.xlabel("Number of samples")
    plt.ylabel("Number of genera (unfiltered)")
    plt.title("Rarefaction curves per compartment")
    plt.legend(title="Compartment")
    plt.tight_layout()

    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def plot_venn_diagram(save_as, figsize=(6,6)):
    compartments = ["Rhizosphere", "Endosphere", "Bulk soil"]
    feature_sets = {}

    for compartment in compartments:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl"
        with open(pkl_file, "rb") as handle:
            ds = pkl.load(handle)

        features = set(ds.taxonomy_counts_df.columns[ds.taxonomy_counts_df.sum(axis=0) > 0])
        feature_sets[compartment] = features

    plt.figure(figsize=figsize)

    colors = config["colors"]["RootCompartment"]
    venn = venn3(
        subsets=(feature_sets["Rhizosphere"],
                 feature_sets["Endosphere"],
                 feature_sets["Bulk soil"]),
        set_labels=compartments,
        set_colors=(colors["Rhizosphere"], colors["Endosphere"], colors["Bulk soil"])
    )

    for subset in ["100", "010", "001"]:
        if venn.get_patch_by_id(subset):
            venn.get_patch_by_id(subset).set_alpha(0.5)

    plt.title("Feature Overlap Between Compartments")
    plt.tight_layout()

    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()

def main():
    # Plot study x study heatmap
    pkl_file = f"../datasets/genus.pkl"
    with open(pkl_file, 'rb') as handle:
        ds = pkl.load(handle)
    plot_jaccard_heatmap(ds,
                         "../plots/jaccard")

    plot_rarefaction_per_compartment("../plots/rarefaction")
    plot_venn_diagram("../plots/venn")

if __name__ == "__main__":
    main()
