import math
import pickle as pkl
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from Bio import Phylo

import yaml

with open('config.yml') as f:
    config = yaml.safe_load(f)

def compute_host_phylogenetic_distances(tree_file):
    tree = Phylo.read(str(tree_file), "newick")
    terminals = tree.get_terminals()
    distance_dict = defaultdict(dict)
    for t1 in terminals:
        name1_space = t1.name.replace("_", " ")
        for t2 in terminals:
            if t1 == t2:
                continue
            name2_space = t2.name.replace("_", " ")
            dist = tree.distance(t1, t2) / 4.0
            distance_dict[name1_space][name2_space] = dist
    return distance_dict


def compute_host_microbiome_similarities(ds,
                                         config,
                                         within_studies=True,
                                         thresh=0.95):
    all_hosts = list(ds.metadata_df["HostSpecific"].unique())
    tax_df = ds.taxonomy_counts_df
    mdf = ds.metadata_df

    similarity_dict = defaultdict(dict)
    cores_dict = {}

    valid_hosts = []
    for host in all_hosts:
        mask = (mdf["HostSpecific"] == host)
        n = int(mask.sum())
        if n < 10:
            continue

        min_count = math.ceil(thresh * n)
        counts = tax_df.loc[mask, :]
        presence_counts = (counts > 0).sum(axis=0)

        core_set = set(presence_counts[presence_counts >= min_count].index)
        cores_dict[host] = core_set
        valid_hosts.append(host)

    for host1 in valid_hosts:
        for host2 in valid_hosts:
            if host1 == host2:
                continue

            pres1_set = cores_dict[host1]
            pres2_set = cores_dict[host2]

            union_size = len(pres1_set | pres2_set)
            if union_size == 0:
                similarity = np.nan
            else:
                intersection_size = len(pres1_set & pres2_set)
                similarity = intersection_size / union_size

            similarity_dict[host1][host2] = similarity

    return similarity_dict, cores_dict


def plot_similarity_vs_phylo(similarity_dict, phylo_distances, title, save_path):
    points_x = []
    points_y = []
    labels = []
    seen_pairs = set()
    
    for h1 in similarity_dict:
        for h2 in similarity_dict[h1]:
            if h1 == h2:
                continue
            pair = tuple(sorted([h1, h2]))
            if pair in seen_pairs:
                continue
            if (h1 not in phylo_distances):
                continue

            seen_pairs.add(pair)

            sim = similarity_dict[h1][h2]
            if np.isnan(sim):
                continue

            if h1 not in phylo_distances:
                continue
            elif h2 not in phylo_distances[h1]:
                continue
            
            phy = phylo_distances[h1][h2]

            points_x.append(sim)
            points_y.append(phy)
            labels.append(f"{pair[0]} | {pair[1]}")

    if not points_x:
        print("Nothing to plot")
        return

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(points_x, points_y, alpha=0.7)
    ax.set_xlabel("Jaccard similarity (microbiome)")
    ax.set_ylabel("Phylogenetic distance")
    ax.set_title(title)
    ax.grid(linestyle='--', alpha=0.4)
    #ax.set_xlim(0, 1)
    plt.tight_layout()

    save_path = Path(save_path)
    fig.savefig(str(save_path.with_suffix(".svg")), format="svg")
    fig.savefig(str(save_path.with_suffix(".png")), format="png", dpi=600)
    plt.close(fig)

def plot_core_sizes(host_dict,
                    save_as):
    hosts = list(host_dict.keys())
    core_sizes = [len(core) for core in host_dict.values()]

    plt.figure(figsize=(0.4*len(hosts), 6))
    plt.bar(hosts, core_sizes, color="steelblue")
    plt.xticks(rotation=90)
    plt.ylabel("Core size (# taxa)")
    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def main():
    dicot_tree = Path("../data/plant_genes/dicots_tree_11.08.2025.nwk")
    mono_tree = Path("../data/plant_genes/monocots_tree_11.08.2025.nwk")
    out_dir = Path("../data/plots/")

    phylo_dicots = compute_host_phylogenetic_distances(dicot_tree)
    phylo_monocots = compute_host_phylogenetic_distances(mono_tree)

    print(phylo_dicots)

    for compartment in ["Rhizosphere", "Endosphere"]:
        pkl_file = f"../data/merged_l6_unfiltered_{compartment.lower().replace(' ', '_')}_batch_corrected.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        similarity, host_dict = compute_host_microbiome_similarities(ds,
                                                                     config,
                                                                     within_studies=False,
                                                                     thresh=0.95)

        plot_similarity_vs_phylo(
            similarity,
            phylo_monocots,
            title=f"Monocots — {compartment}",
            save_path=out_dir / f"monocots_{compartment.replace(' ','_')}_sim_vs_phylo"
        )

        plot_similarity_vs_phylo(
            similarity,
            phylo_dicots,
            title=f"Dicots — {compartment}",
            save_path=out_dir / f"dicots_{compartment.replace(' ','_')}_sim_vs_phylo"
        )

        plot_core_sizes(host_dict,
                        save_as=f"../data/plots/bar_core_sizes_{compartment.replace(' ','_')}")


if __name__ == "__main__":
    main()
