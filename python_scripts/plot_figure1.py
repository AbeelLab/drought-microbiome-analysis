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
from utils import trim_taxonomy

with open('config.yml') as f:
    config = yaml.safe_load(f)

def extract_level(tax_str):
        for part in tax_str.split(';'):
            if part.startswith("p__"):
                return part.replace("p__", '')
        return 'Unassigned'

def plot_taxonomy(ds,
                  save_as,
                  figsize=(10,4),
                  ordered_hosts=["PACMAD clade", "BOP clade", "Fabids", "Asterids", "Malvids"],
                  level="p__"):
    from copy import deepcopy
    import pandas as pd
    import matplotlib.pyplot as plt

    ds_copy = deepcopy(ds)
    normalized = ds_copy.get_normalized_features()
    colors = config["colors"]["Phyla"]

    phyla_order = list(colors.keys())  # Use all phyla in config as plotting order

    results = []

    for host in ordered_hosts:
        host_samples = ds.metadata_df[ds.metadata_df["Host"] == host].index
        normalized_host_df = normalized.loc[host_samples]

        phylum_labels = [extract_level(c) for c in normalized_host_df.columns]
        normalized_host_df.columns = phylum_labels
        phylum_abund = normalized_host_df.groupby(normalized_host_df.columns, axis=1).sum()

        phylum_mean = phylum_abund.mean(axis=0)

        # Iterate through all phyla in config
        for phylum in phyla_order:
            if phylum in phylum_mean:
                results.append({"Host": host, "Phylum": phylum, "Abundance": phylum_mean[phylum]})
        
        # Sum everything else outside config as "Other"
        other_sum = phylum_mean.drop(phyla_order, errors='ignore').sum()
        if other_sum > 0:
            results.append({"Host": host, "Phylum": "Other", "Abundance": other_sum})

    plot_df = pd.DataFrame(results)

    plt.figure(figsize=figsize)
    left = {h: 0 for h in ordered_hosts}
    legend_added = set()

    for host in ordered_hosts:
        host_df = plot_df[plot_df["Host"] == host].copy()
        # Sort according to config order, 'Other' last
        host_df["sort_order"] = host_df["Phylum"].apply(lambda x: phyla_order.index(x) if x in phyla_order else float('inf'))
        host_df = host_df.sort_values(by="sort_order")
        for _, row in host_df.iterrows():
            label = row["Phylum"] if row["Phylum"] not in legend_added else None
            plt.barh(row["Host"],
                     row["Abundance"],
                     left=left[host],
                     color=colors.get(row["Phylum"], "gray"),
                     label=label)
            left[host] += row["Abundance"]
            legend_added.add(row["Phylum"])

    plt.xlabel("Relative abundance (%)")
    plt.ylabel("Host")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", title="Phylum")
    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()
    

def plot_feature_distribution(ds, save_as, figsize=(8, 5), bins=100):
    sample_sums = ds.taxonomy_counts_df.sum(axis=1)

    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(sample_sums, bins=bins, color="darkorange", edgecolor="black", alpha=0.7)
    ax.set_xlabel("Total Feature Counts per Sample")
    ax.set_ylabel("Number of Samples")
    ax.set_title("Distribution of Total Counts per Sample")
    ax.grid(axis="y", linestyle="--", alpha=0.5)

    plt.tight_layout()
    fig.savefig(f"{save_as}.svg", format="svg")
    fig.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close(fig)

    return sample_sums

def plot_core_thresholds(dataset,
                         save_as="../data/plots/core_thresholds",
                         step=0.05,
                         figsize=(6, 7),
                         color=None,
                         figsize_base_specific=(8, 4)):
    thresholds = np.arange(step, 1.0, step)
    counts = [len(dataset.get_core(t)[0]) for t in thresholds]

    fig, ax = plt.subplots(figsize=figsize)
    ax.plot(thresholds * 100, counts, marker='o', linewidth=2, color=color)
    ax.set_xlabel('Prevalence threshold (%)')
    ax.set_ylabel('Number of core taxa')
    ax.set_title('Core taxa vs. prevalence threshold')
    ax.set_xticks(thresholds * 100)
    ax.set_ylim(0, 650)
    ax.grid(axis='y', linestyle='--', alpha=0.5)

    for t, c in zip(thresholds, counts):
        if t >= 0.75:
            ax.text(t*100, c, str(c), va='bottom', ha='center', fontsize=9)

    plt.tight_layout()

    if save_as:
        fig.savefig(f"{save_as}.svg", format='svg')
        fig.savefig(f"{save_as}.png", format='png', dpi=600)
    plt.close(fig)

    taxa_above_90 = dataset.get_core(0.9)[0]
    return taxa_above_90
    

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

def plot_rarefaction_per_compartment(save_as, n_permutations=100, figsize=(8, 5)):
    plt.figure(figsize=figsize)

    np.random.seed(42)
    for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
        pkl_file = f"../data/merged_l6.pkl"
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
        pkl_file = f"../data/merged_l6_unfiltered_{compartment.lower().replace(' ', '_')}.pkl"
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


def plot_core_taxa(ds,
                   core,
                   save_as,
                   level="p__"):

    colors = config["colors"]["Phyla"]
    norm = ds.get_normalized_features()

    def extract_level(tax_str):
        for part in tax_str.split(';'):
            if part.startswith(level):
                return part.replace(level, '')
        return 'Unassigned'

    records = []
    for taxon in core:
        abundances = norm[taxon]
        nonzero = abundances[abundances > 0]
        phylum = extract_level(taxon)
        color = colors.get(phylum, "#E8E9EB")
        for val in nonzero:
            records.append({"Taxon": trim_taxonomy(taxon), "Abundance": val,
                            "Phylum": phylum, "Color": color})

    df = pd.DataFrame(records)
    taxa_order = df.groupby("Taxon")["Abundance"].mean().sort_values(ascending=False).index

    # Width proportional to number of core taxa
    figsize = (len(core) * 0.8, 10)

    plt.figure(figsize=figsize)
    palette = {p: colors.get(p, "#E8E9EB") for p in df["Phylum"].unique()}

    ax = sns.boxplot(
        data=df,
        x="Taxon",
        y="Abundance",
        order=taxa_order,
        showcaps=True,
        showbox=True,
        showfliers=True,
        hue="Phylum",
        whiskerprops={'linewidth': 2},
        medianprops={'linewidth': 2},
        fill=False,
        palette=palette
    )

    ax.set_ylim(-1, 38)  # y-axis limits

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title("Core taxa non-zero abundances")
    plt.ylabel("Relative abundance (%)")
    plt.xlabel("Taxon")
    plt.tight_layout()

    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def main():
    # Plot study x study heatmap
    pkl_file = f"../data/merged_l6_filtered.pkl"
    with open(pkl_file, 'rb') as handle:
        ds = pkl.load(handle)
    plot_jaccard_heatmap(ds,
                         "../data/plots/jaccard")

    plot_rarefaction_per_compartment("../data/plots/rarefaction")
    plot_venn_diagram("../data/plots/venn")
    

    for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
        pkl_file = f"../data/merged_l6_filtered_{compartment.lower().replace(' ', '_')}.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        # Plot dataset statistics
        core = plot_core_thresholds(ds,
                                    save_as=f"../data/plots/core_thresholds_{compartment}",
                                    color=config["colors"]["RootCompartment"][compartment])
        plot_core_taxa(ds,
                       core,
                       save_as=f"../data/plots/core_box_plot_thresholds_{compartment}")
        
        plot_feature_distribution(ds,
                                  f"../data/plots/feature_distribution_{compartment}")
        
        if compartment != "Bulk soil":
            plot_taxonomy(ds,
                          save_as=f"../data/plots/taxonomy_{compartment}")
        else:
            plot_taxonomy(ds,
                          figsize=(10, 2),
                          ordered_hosts=["Soil"],
                          save_as=f"../data/plots/taxonomy_{compartment}")


if __name__ == "__main__":
    main()

