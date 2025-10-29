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


def collect_core_data(ds, core, level="p__"):
    colors = config["colors"]["Phyla"]
    norm = ds.get_normalized_features()

    records = []
    for taxon in core:
        abundances = norm[taxon]
        nonzero = abundances[abundances > 0]
        phylum = extract_level(taxon)
        color = colors.get(phylum, "#E8E9EB")
        for val in nonzero:
            records.append({
                "Taxon": trim_taxonomy(taxon),
                "Abundance": val,
                "Phylum": phylum,
                "Color": color
            })
    return pd.DataFrame(records)

def plot_core_taxa_combined(compartments, save_as="../plots/core_boxplot_combined.png"):
    dfs = []
    taxa_order = []

    for compartment in compartments:
        with open(f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl", "rb") as handle:
            ds = pkl.load(handle)

        core = ds.get_core(0.9)[0]
        df = collect_core_data(ds, core)
        df["Compartment"] = compartment
        dfs.append(df)

        for taxon in df["Taxon"].unique():
            if taxon not in taxa_order:
                taxa_order.append(taxon)

    # Combine all dataframes
    full_df = pd.concat(dfs, ignore_index=True)

    # Palette based on config colors
    colors = config["colors"]["Phyla"]
    palette = {p: colors.get(p, "#E8E9EB") for p in full_df["Phylum"].unique()}

    # Consistent box width
    box_width_inch = 0.5
    base_padding_inch = 2
    unique_taxa = len(taxa_order)
    base_fig_width = unique_taxa * box_width_inch + base_padding_inch

    # Define subplot heights based on expected y-limits
    heights = []
    for c in compartments:
        if c.lower() == "endosphere":
            heights.append(1.5)  # taller panel
        else:
            heights.append(1.0)

    # Normalize heights
    total_height = sum(heights)
    fig_height = 3.5 * total_height

    fig, axes = plt.subplots(
        nrows=len(compartments),
        ncols=1,
        figsize=(base_fig_width, fig_height),
        sharex=True,
        sharey=False,
        gridspec_kw={'height_ratios': heights}
    )

    if len(compartments) == 1:
        axes = [axes]

    for ax, compartment in zip(axes, compartments):
        df = full_df[full_df["Compartment"] == compartment]
        sns.boxplot(
            data=df,
            x="Taxon",
            y="Abundance",
            order=taxa_order,
            hue="Phylum",
            showcaps=True,
            showfliers=True,
            whiskerprops={'linewidth': 2},
            medianprops={'linewidth': 2},
            fill=False,
            palette=palette,
            width=0.6,
            ax=ax
        )

        ax.set_title(compartment)

        # Conditional y-limits
        if compartment.lower() == "endosphere":
            ax.set_ylim(-1, 65)
        else:
            ax.set_ylim(-1, 30)

        ax.set_ylabel("Relative abundance (%)")
        ax.legend([], [], frameon=False)

    axes[-1].set_xticklabels(axes[-1].get_xticklabels(), rotation=90)
    axes[-1].set_xlabel("Taxon")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", title="Phylum")

    plt.tight_layout(rect=[0, 0, 0.98, 1])
    fig.savefig(save_as, dpi=600, bbox_inches="tight")
    fig.savefig(save_as.replace(".png", ".svg"), format="svg")
    plt.close(fig)

    

def main():
    compartments = ["Endosphere", "Rhizosphere", "Bulk soil"]
    plot_core_taxa_combined(compartments)
    
    for compartment in compartments:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        # Plot dataset statistics
        core = plot_core_thresholds(ds,
                                    save_as=f"../plots/core_thresholds_{compartment}",
                                    color=config["colors"]["RootCompartment"][compartment])

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
