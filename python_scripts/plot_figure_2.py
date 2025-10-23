import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pickle as pkl
import seaborn as sns
import yaml

from copy import deepcopy
from matplotlib.patches import Patch
from utils import trim_taxonomy

with open('config.yml') as f:
    config = yaml.safe_load(f)

def extract_level(tax_str, level="p__"):
        for part in tax_str.split(';'):
            if part.startswith(level):
                return part.replace(level, '')
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

        phylum_mean = normalized_host_df.mean(axis=0)

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


def plot_core_taxa(ds,
                   core,
                   save_as,
                   level="p__"):

    colors = config["colors"]["Phyla"]
    norm = ds.get_normalized_features()

    print(core)
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

    ax.set_ylim(-1, 38)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.title("Core taxa non-zero abundances")
    plt.ylabel("Relative abundance (%)")
    plt.xlabel("Taxon")
    plt.tight_layout()

    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def main():
    for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        # Plot dataset statistics
        core = plot_core_thresholds(ds,
                                    save_as=f"../plots/core_thresholds_{compartment}",
                                    color=config["colors"]["RootCompartment"][compartment])
        plot_core_taxa(ds,
                       core,
                       save_as=f"../plots/core_box_plot_thresholds_{compartment}")

        # For taxonomic composition, use phylum-level
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/phylum.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)
        
        if compartment != "Bulk soil":
            plot_taxonomy(ds,
                          save_as=f"../plots/taxonomy_{compartment}")
        else:
            plot_taxonomy(ds,
                          figsize=(10, 2),
                          ordered_hosts=["Soil"],
                          save_as=f"../plots/taxonomy_{compartment}")


if __name__ == "__main__":
    main()
