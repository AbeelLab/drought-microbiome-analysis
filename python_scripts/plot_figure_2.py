import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import pickle as pkl
import seaborn as sns
import os
import yaml

from analysis import run_mmuphin_diff_abundance
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


def plot_core_taxa_combined(compartments, save_as="../plots/core_scatter_combined.png"):
    dfs = []

    # Gather data for all compartments
    for compartment in compartments:
        with open(f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl", "rb") as handle:
            ds = pkl.load(handle)

        data = []
        for taxon in ds.taxonomy_counts_df.columns:
            if "uncultured" in taxon.lower():
                continue

            prevalence = ds.get_prevalence(taxon)
            if prevalence < 80:
                continue

            mean_abundance = ds.get_normalized_features()[taxon][ds.taxonomy_counts_df[taxon] > 0].mean()
            phylum = extract_level(taxon)

            data.append({
                "Taxon": taxon,
                "Prevalence": prevalence,
                "MeanAbundance": mean_abundance,
                "Phylum": phylum,
                "Compartment": compartment
            })

        if data:
            dfs.append(pd.DataFrame(data))

    if not dfs:
        print("[WARNING] No taxa passed the filtering criteria")
        return

    full_df = pd.concat(dfs, ignore_index=True)

    # Count in how many compartments each taxon passes the 80% prevalence
    taxon_compartment_counts = (
        full_df.groupby("Taxon")["Compartment"]
        .nunique()
        .to_dict()
    )

    # Map number of compartments to point size
    def map_size(n):
        if n == 1:
            return 60   # small
        elif n == 2:
            return 180  # medium
        else:
            return 360  # big

    full_df["Size"] = full_df["Taxon"].map(lambda t: map_size(taxon_compartment_counts.get(t, 1)))

    # Map phyla to colors
    colors = config["colors"]["Phyla"]
    full_df["Color"] = full_df["Phylum"].map(lambda p: colors.get(p, "#E8E9EB"))

    # Create figure
    fig, axes = plt.subplots(
        nrows=len(compartments),
        ncols=1,
        figsize=(7.5, 4*len(compartments)),
        sharex=True
    )
    if len(compartments) == 1:
        axes = [axes]

    # Scatter plot per compartment
    for ax, compartment in zip(axes, compartments):
        df_c = full_df[full_df["Compartment"] == compartment]

        ax.scatter(
            df_c["MeanAbundance"],
            df_c["Prevalence"],
            c=df_c["Color"],
            s=df_c["Size"],
            edgecolor="k"
        )

        # Label each point
        for _, row in df_c.iterrows():
            ax.text(
                row["MeanAbundance"] + 0.05,
                row["Prevalence"] + 0.05,
                extract_level(row["Taxon"], "g__"),
                fontsize=12
            )

        ax.set_title(compartment, fontsize=16)
        ax.set_ylabel("Prevalence (%)", fontsize=16)
        ax.set_xlabel("Mean non-zero abundance (%)", fontsize=16)
        ax.set_xlim(left=0)
        ax.set_ylim(75, 100)

        # Increase tick label sizes
        ax.tick_params(axis='both', which='major', labelsize=14)

    plt.tight_layout()
    fig.savefig(save_as, dpi=600, bbox_inches="tight")
    fig.savefig(save_as.replace(".png", ".svg"), format="svg")
    plt.close(fig)



def main():
    compartments = ["Endosphere", "Rhizosphere", "Bulk soil"]
    plot_core_taxa_combined(compartments)

    # DA analysis between phyla in different compartments
    pkl_file = f"../datasets/phylum_batch_corrected.pkl"
    with open(pkl_file, 'rb') as handle:
        batch_corrected_phyla = pkl.load(handle)

    for compartment in compartments:
        ds_copy = deepcopy(batch_corrected_phyla)
        ds_copy.filter_samples_based_on_condition(lambda df: df["RootCompartment"] != compartment)
        ds_copy.remove_zero_features()

        results, ref = run_mmuphin_diff_abundance(ds_copy, exposure="RootCompartment")

        df_results = pd.DataFrame.from_dict(results, orient="index")[["logFC", "adj.P.Val"]]
        df_results = df_results.sort_values("adj.P.Val", ascending=True)

        other_levels = [x for x in compartments if x != compartment]
        comparison_label = f"{' vs '.join(other_levels)} (reference: {ref})"
        filename = f"DA_{'_vs_'.join(other_levels).replace(' ', '_')}_ref_{ref.replace(' ', '_')}.csv"
        filepath = os.path.join("../raw_data", filename)

        with open(filepath, "w") as f:
            with open(filepath, "w") as f:
                f.write(f"# {comparison_label}\n")
            df_results.to_csv(filepath, mode="a", header=True)

        print(f"[INFO] Saved differential abundance results to {filepath}")
    
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

        # Visualize: no batch correction
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
