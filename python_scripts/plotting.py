import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import seaborn as sns
import os
import numpy as np
import re
import textwrap
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist, squareform
import math
from matplotlib.patches import Patch
from scipy.stats import spearmanr
import networkx as nx
import yaml
from scipy.stats import spearmanr
from matplotlib.patches import Rectangle
from copy import deepcopy

with open('config.yml') as f:
    config = yaml.safe_load(f)
color_map = config["colors"]

def trim_taxonomy(tax_str):
    parts = tax_str.split(";")
    for part in reversed(parts):
        part = part.strip()
        if part not in ["", "__"]:
            return part
    return ""


def plot_corr(phy_distances, similarities, dot_color, alpha, save_as):
    x, y, seen = [], [], set()
    for a, inner in similarities.items():
        for b, sim in inner.items():
            if np.isnan(sim):
                continue
            pair = tuple(sorted((a, b)))
            if pair in seen:
                continue
            seen.add(pair)
            phy = phy_distances.get(a, {}).get(b) or phy_distances.get(b, {}).get(a)
            if phy is None:
                continue
            x.append(sim)
            y.append(phy)

    rho, pval = spearmanr(x, y)
    n_points = len(x)

    fig, ax = plt.subplots()
    ax.scatter(x, y,
               facecolors=dot_color,
               edgecolors='black',
               linewidths=0.5,
               alpha=alpha)
    
    ax.set_xlabel("# common genera")
    ax.set_ylabel("Phylogenetic distance")
    ax.set_title(f"{n_points} pairs — Spearman ρ={rho:.2f}, p={pval:.2g}")
    ax.set_xlim(0, 500)
    ax.set_ylim(0, 0.15)
    plt.tight_layout()

    for ext in ('svg', 'png'):
        fig.savefig(f"{save_as}.{ext}", dpi=600)
    plt.close(fig)

def plot_pcoa(dataset,
              save_as):
    X = dataset.get_counts_features()
    metadata = dataset.get_metadata()

    bc_distances = pdist(X.values, metric='braycurtis')
    distance_matrix = DistanceMatrix(squareform(bc_distances), ids=X.index)

    pcoa_object = pcoa(distance_matrix)
    pcoa_result = pcoa_object.samples[["PC1", "PC2"]]
    pcoa_df = pd.DataFrame(pcoa_result, columns=["PC1", "PC2"], index=X.index)

    # Join metadata (other features)
    pcoa_df = pd.concat([pcoa_df, metadata], axis=1)

    def save_pcoa_plot(color_by):
        plt.figure(figsize=(8, 6))

        sns.scatterplot(data=pcoa_df,
                        x="PC1", y="PC2",
                        hue=color_by,
                        palette=color_map[color_by],
                        alpha=0.75)

        plt.title(f"PCoA: PC1 vs PC2 (Bray-Curtis)\nColored by: {color_by}")
        plt.legend(title=color_by, bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.savefig(save_as + "_" + color_by + ".svg", bbox_inches="tight")
        plt.savefig(save_as + "_" + color_by + ".png", dpi=600, bbox_inches="tight")
        print(f"[INFO] PCoA plot saved to {save_as}_{color_by}")
        plt.close()

    for category in ["StudyName", "Host"]:
        save_pcoa_plot(category)

        
def plot_heatmap(ds,
                 ordered_ticks):
    to_compare = []
    for host, study in ordered_ticks:
        mask = (ds.metadata_df["Host"] == host) & (ds.metadata_df["StudyName"] == study)
        sub_df = ds.taxonomy_counts_df.loc[mask]
        nonzero_cols = sub_df.columns[(sub_df.sum(axis=0) > 0).values]
        to_compare.append(set(nonzero_cols))

    n = len(to_compare)
    overlap = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            inter = len(to_compare[i].intersection(to_compare[j]))
            denom = len(to_compare[i]) + len(to_compare[j])
            overlap[i, j] = (2 * inter / denom) if denom > 0 else 0

    print(overlap)
    fig, ax = plt.subplots()
    im = ax.imshow(overlap,
                   cmap='Greys',
                   vmin=0,
                   vmax=1,
                   aspect='auto')

    labels = [f"{host}+{study}" for host, study in ordered_ticks]
    ax.set_xticks(np.arange(n))
    ax.set_yticks(np.arange(n))
    ax.set_xticklabels(labels, rotation=90)
    ax.set_yticklabels(labels)

    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Shared Features Count')

    fig.savefig("../data/plots/" + ds.dataset_name + "_heatmap.png", dpi=300, bbox_inches='tight')
    fig.savefig("../data/plots/" + ds.dataset_name + "_heatmap.svg", bbox_inches='tight')
    plt.close()

def plot_taxonomy(ds,
                  ordered_hosts,
                  colors,
                  save_as,
                  level="p__"):
    ds_copy = deepcopy(ds)
    normalized = ds_copy.get_normalized_features()

    def extract_level(tax_str):
        for part in tax_str.split(';'):
            if part.startswith(level):
                return part.replace(level, '')
        return 'Unassigned'

    # Extract taxonomic levels
    phylum_labels = [extract_level(c) for c in normalized.columns]
    normalized.columns = phylum_labels
    phylum_abund = normalized.groupby(normalized.columns, axis=1).sum()

    # Top phyla + "Other"
    mean_abund = phylum_abund.mean(axis=0)
    top6 = mean_abund.nlargest(6).index.tolist()
    phylum_top = phylum_abund[top6].copy()
    phylum_top['Other'] = 100 - phylum_top.sum(axis=1)

    # Grouping by Host | StudyName
    md = ds_copy.metadata_df.copy()
    md['Group'] = md['Host'].astype(str) + ' | ' + md['StudyName'].astype(str)
    df = phylum_top.join(md['Group'])
    group_means = df.groupby('Group')[top6 + ['Other']].mean()

    # Reindex to ordered_hosts and fill missing
    group_means = group_means.reindex(ordered_hosts).fillna(0)

    # Prepare spacing and centered dashed lines
    hosts = [label.split(' | ')[0] for label in group_means.index]
    y_positions = []
    y_tick_labels = []
    current_y = 0
    host_line_positions = []

    for i, label in enumerate(group_means.index):
        if i > 0 and hosts[i] != hosts[i - 1]:
            current_y += 1  # Add space between hosts
            host_line_positions.append(current_y - 0.5)  # Center of the added space
        y_positions.append(current_y)
        y_tick_labels.append(label)
        current_y += 1

    height = max(4, current_y * 0.3)
    fig, ax = plt.subplots(figsize=(10, height))

    # Draw stacked bars
    cumulative = np.zeros(len(y_positions))
    for phylum in group_means.columns:
        vals = group_means[phylum].values
        color = colors.get(phylum, 'gray')
        ax.barh(y_positions, vals, left=cumulative, label=phylum, color=color)
        cumulative += vals

    # Draw centered dashed lines between hosts
    for pos in host_line_positions:
        ax.axhline(pos, color='black', linestyle='--', linewidth=1.5, alpha=0.8)

    # Axes formatting
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_tick_labels)
    ax.set_xlabel('Mean Relative Abundance (%)')
    ax.set_ylabel('Host | StudyName')
    ax.set_xlim(0, 100)
    ax.set_ylim(-0.5, current_y - 0.5)

    # Full legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=colors[p], label=p) for p in colors]
    ax.legend(handles=handles, title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    fig.savefig(f"{save_as}.svg", format='svg')
    fig.savefig(f"{save_as}.png", format='png', dpi=600)
    plt.close(fig)



def plot_num_samples(ds,
                     ordered_hosts,
                     save_as):
    metadata = ds.metadata_df.copy()
    metadata['Host_Study'] = metadata['Host'] + ' | ' + metadata['StudyName']
    metadata = metadata[metadata['Host_Study'].isin(ordered_hosts)]
    metadata['Host_Study'] = pd.Categorical(metadata['Host_Study'], categories=ordered_hosts, ordered=True)

    counts = metadata.groupby(['Host_Study', 'Treatment']).size().unstack(fill_value=0)
    counts = counts.loc[ordered_hosts]

    # Prepare bar data
    drought = counts.get('Drought', pd.Series(0, index=counts.index))
    control = counts.get('Control', pd.Series(0, index=counts.index))

    # Build y-axis with extra space between host groups
    hosts = [label.split(' | ')[0] for label in counts.index]
    y_positions = []
    y_tick_labels = []
    current_y = 0
    host_boundaries = []

    for i, label in enumerate(counts.index):
        if i > 0 and hosts[i] != hosts[i - 1]:
            current_y += 1  # Extra space between different hosts
            host_boundaries.append(current_y - 0.5)
        y_positions.append(current_y)
        y_tick_labels.append(label)
        current_y += 1

    height = max(4, current_y * 0.4)
    fig, ax = plt.subplots(figsize=(10, height))

    # Plot bars with spacing
    ax.barh(y_positions, drought.values, edgecolor="white", color="#fd9a2e", label="Drought")
    ax.barh(y_positions, control.values, edgecolor="white", left=drought.values, color="#3098fe", label="Control")

    # Draw thick dashed lines between host groups
    for pos in host_boundaries:
        ax.axhline(pos, color='black', linestyle='--', linewidth=1.5, alpha=0.8)

    # Format axes
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_tick_labels)
    ax.set_xlabel('Number of Samples')
    ax.set_ylabel('Host | StudyName')
    ax.set_title('Number of Samples per Host and Study by Treatment')
    ax.set_ylim(-0.5, current_y - 0.5)
    ax.legend()

    plt.tight_layout()
    fig.savefig(f"{save_as}.svg", format='svg')
    fig.savefig(f"{save_as}.png", format='png', dpi=600)
    plt.close(fig)


def plot_core_taxa(ds,
                   ordered_hosts,
                   core_taxa,
                   save_as):
    metadata = ds.metadata_df
    taxonomy = ds.taxonomy_counts_df

    metadata['Host_Study'] = metadata['Host'] + ' | ' + metadata['StudyName']

    core_counts = []
    non_core_counts = []

    for host_study in ordered_hosts:
        sample_ids = metadata[metadata['Host_Study'] == host_study].index
        subset = taxonomy[sample_ids]
        taxa_present = subset[(subset > 0).any(axis=1)].index
        n_core = len(set(taxa_present) & set(core_taxa))
        n_total = len(taxa_present)
        n_non_core = n_total - n_core
        core_counts.append(n_core)
        non_core_counts.append(n_non_core)

    y = range(len(ordered_hosts))
    plt.figure(figsize=(10, len(ordered_hosts) * 0.4))

    plt.barh(y, core_counts, edgecolor="white", color='black', label='Core taxa')
    plt.barh(y, non_core_counts, left=core_counts, edgecolor="white", color='gray', label='Non-core taxa')

    plt.yticks(y, ordered_hosts)
    plt.xlabel('Number of Taxa')
    plt.title('Core vs Non-Core Taxa per Host | StudyName')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format='svg')
    plt.savefig(f"{save_as}.png", format='png', dpi=600)
    plt.close()
    


# Level doesn't matter actually and the plot should be the same for all levels
# Because it is a sample plot
# Set phylum level as default
def plot_bubbles(merged_file, plotting_dir, level=2):
    region_colors = {
        'California':   '#01665e',
        'Michigan':     '#35978f',
        'Saskatchewan': '#80cdc1',
        'Ontario':      '#c7eae5',
    }
    soil_colors = {
        'Clay':  '#dfc27d',
        'Sand':  '#bf812d',
    }

    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    loc_counts  = df[df['Location'] != 'Greenhouse'].groupby('Location').size()
    soil_counts = df[df['Location'] == 'Greenhouse'].groupby('SoilType').size()

    all_items = []
    max_count = max(loc_counts.max(), soil_counts.max())

    for i, (loc, count) in enumerate(loc_counts.items()):
        all_items.append({
            'label': loc,
            'count': count,
            'color': region_colors.get(loc, '#000000'),
            'group': 'Region'
        })

    for i, (soil, count) in enumerate(soil_counts.items()):
        all_items.append({
            'label': soil,
            'count': count,
            'color': soil_colors.get(soil, '#000000'),
            'group': 'Soil'
        })

    # Plot
    fig, ax = plt.subplots(figsize=(len(all_items)*1.5, 4), dpi=150)

    x = range(len(all_items))
    y = [0] * len(all_items)
    sizes = [(item['count'] * 10) for item in all_items]  # scale bubble area

    for i, item in enumerate(all_items):
        ax.scatter(x[i], y[i], s=sizes[i], color=item['color'], edgecolors='black', alpha=0.8)
        ax.text(x[i], 0.1, f"{item['label']} ({item['count']})",
                ha='center', va='bottom', fontsize=10)

    ax.axis('off')
    ax.set_xlim(-1, len(all_items))
    ax.set_ylim(-1, 1)

    # Save
    out_file = os.path.join(plotting_dir, "sample_bubbles_scaled.svg")
    plt.savefig(out_file, dpi=600, transparent=True, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Bubble PNG saved to {out_file}")

def plot_lfc_diff_abundance(diff_abundance_results,
                            level,
                            config):
    plotting_dir = config["plotting_dir"]
    print("Microbial signature size:", len(diff_abundance_results))

    # Prepare results with phylum and only cultured taxa
    records = []
    for full_tax, metrics in diff_abundance_results.items():
        # Extract phylum (string after 'p__' and before next ';')
        match = re.search(r'p__([^;]+)', full_tax)
        phylum = match.group(1) if match else 'Unknown'
        name = full_tax.split(';')[-1]  # last taxonomic rank
        name = name.replace('s__', '').replace('g__', '').replace('f__', '') \
                   .replace('o__', '').replace('c__', '').replace('p__', '')
        if name.isalpha():  # cultured taxa filter
            records.append({
                'taxon': name,
                'logFC': metrics['logFC'],
                'phylum': phylum
            })

    df = pd.DataFrame(records)

    # Sort and select top/bottom 10 by logFC
    df_sorted = df.sort_values('logFC', ascending=False)
    top10 = df_sorted.head(10)
    bottom10 = df_sorted.tail(10)
    plot_df = pd.concat([top10, bottom10])

    # Group by phylum: sort within each phylum by logFC descending
    plot_df['phylum'] = plot_df['phylum'].astype(str)
    plot_df = plot_df.sort_values(['phylum', 'logFC'], ascending=[True, False])

    # Compute bar positions with extra gap between positive and negative groups
    # Identify index split between positive and negative values after grouping
    positive_df = plot_df[plot_df['logFC'] >= 0]
    negative_df = plot_df[plot_df['logFC'] < 0]

    n_pos = len(positive_df)
    n_neg = len(negative_df)
    gap = 3.5  # multiplier for extra spacing between the two sign groups
    pos_positions = np.arange(n_pos)
    neg_positions = np.arange(n_pos + gap, n_pos + gap + n_neg)

    # Combine positions preserving order: positives first, then negatives
    positions = np.concatenate([pos_positions, neg_positions])
    combined_df = pd.concat([positive_df, negative_df])

    # Create color map for phyla
    phyla = combined_df['phylum'].unique()
    cmap = plt.get_cmap('tab20')
    color_map = {ph: cmap(i % cmap.N) for i, ph in enumerate(phyla)}
    bar_colors = [color_map[ph] for ph in combined_df['phylum']]

    # Plot
    plt.figure(figsize=(14, 8))
    bars = plt.bar(positions, combined_df['logFC'], color=bar_colors)
    plt.axhline(0, color='gray', linewidth=1)

    # Legend: one handle per phylum
    handles = [plt.Line2D([0], [0], color=color_map[ph], lw=6) for ph in phyla]
    plt.legend(handles, phyla, title='Phylum', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Adjust ticks and labels
    plt.xticks(positions, combined_df['taxon'], rotation=45, ha='right', fontsize=16)
    plt.yticks(fontsize=16)

    plt.title(f"Top/Bottom Log Fold Changes (Level {level})", fontsize=16)
    plt.xlabel('Taxon', fontsize=14)
    plt.ylabel('Log2 Fold Change', fontsize=14)
    plt.tight_layout()

    # Save
    out_file = os.path.join(plotting_dir, f"lfc_l{level}.svg")
    plt.savefig(out_file, dpi=600, transparent=True, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Log fold change plot saved to {out_file}")

    return combined_df['taxon']


def plot_diff_abundance_comparison(plotted_labels,
                                   diff_abundance_results,
                                   diff_abundance_results_inoculum):
    # Collect all taxa short names from main and studies
    def short_name(full_tax):
        name = full_tax.split(';')[-1]
        return re.sub(r'^[a-z]__', '', name)

    main_taxa = {short_name(ft): ft for ft in diff_abundance_results}
    inoc_taxa = set()
    for study_res in diff_abundance_results_inoculum.values():
        for ft in study_res:
            inoc_taxa.add(short_name(ft))

    # Intersection: taxa present in main and at least one study
    common = [tax for tax in main_taxa if tax in inoc_taxa]

    # Order by main logFC descending
    main_logfc = np.array([diff_abundance_results[main_taxa[t]]['logFC'] for t in common])
    order = np.argsort(-main_logfc)
    ordered_taxa = [common[i] for i in order]
    ordered_logfc = main_logfc[order]

    # Determine positions with gap logic
    signs = ordered_logfc >= 0
    n_pos = signs.sum()
    n_total = len(ordered_taxa)
    gap = 3.5
    pos_positions = np.arange(n_pos)
    neg_positions = np.arange(n_pos + gap, n_pos + gap + (n_total - n_pos))
    positions = np.empty(n_total)
    positions[signs] = pos_positions
    positions[~signs] = neg_positions

    # Plot bars
    plt.figure(figsize=(14, 8))
    plt.bar(positions, ordered_logfc, color='lightgray')
    plt.axhline(0, color='gray', linewidth=1)

    # Overlay inoculum points
    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
    for i, (study, results) in enumerate(diff_abundance_results_inoculum.items()):
        xs, ys = [], []
        for idx, tax in enumerate(ordered_taxa):
            full = main_taxa[tax]
            if full in results:
                xs.append(positions[idx])
                ys.append(results[full]['logFC'])
        if xs:
            plt.scatter(xs, ys, marker=markers[i % len(markers)], s=100, label=study)

    plt.legend(title='Inoculum Study', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(positions, ordered_taxa, rotation=45, ha='right', fontsize=16)
    plt.yticks(fontsize=16)
    plt.title('Comparison of Log2 Fold Changes Across Studies', fontsize=16)
    plt.xlabel('Taxon', fontsize=14)
    plt.ylabel('Log2 Fold Change', fontsize=14)
    plt.tight_layout()

    out_file = 'diff_abundance_comparison.svg'
    plt.savefig(out_file, dpi=600, transparent=True, bbox_inches='tight')
    plt.close()
    print(f"[INFO] Comparison plot saved to {out_file}")

    return ordered_taxa, positions

# inputs are two dictionaries with taxonomy names as keys
# and values as dictionries with logFC and p-value
# taxa: {'logFC': ..., 'adj.P.Val': ...}
def make_volcano_plot(diff_abundance_results_inoculum,
                      diff_abundance_results,
                      save_as,
                      total_num_features=None):
    df = pd.DataFrame.from_dict(diff_abundance_results_inoculum, orient='index')
    df.dropna(subset=['logFC', 'adj.P.Val'], inplace=True)
    df['-log10(padj)'] = -np.log10(df['adj.P.Val'])

    df['color'] = 'gray'
    df['size'] = 40
    for taxon, vals in diff_abundance_results.items():
        if taxon in df.index and (lf := vals.get('logFC')) is not None:
            df.at[taxon, 'color'] = '#ff863d' if lf > 0 else 'black'
            df.at[taxon, 'size'] = 60

    # Main plot with colored points
    plt.figure(figsize=(6, 8))
    sns.scatterplot(data=df[df.color=='gray'], x='logFC', y='-log10(padj)',
                    color='gray', alpha=0.5, s=40, edgecolor=None, legend=False)
    sns.scatterplot(data=df[df.color!='gray'], x='logFC', y='-log10(padj)',
                    hue='color', palette={'black':'black','#ff863d':'#ff863d'},
                    size='size', sizes=(60,60), legend=False, edgecolor=None)

    # Stats in title if requested
    if total_num_features is not None:
        sig = df[(df['adj.P.Val'] < 0.05) & (df['logFC'].abs() > 2)]
        n_sig = len(sig)
        n_col = sig['color'].ne('gray').sum()
        pct_sig = n_sig / total_num_features if total_num_features else 0
        pct_col_sig = n_col / n_sig if n_sig else 0
        title = (
            f"Volcano Plot | Total: {total_num_features} | "
            f"Sig: {n_sig} ({pct_sig:.1%}) | "
            f"Colored (significant): {n_col} ({pct_col_sig:.1%})"
        )
    else:
        title = "Volcano Plot"

    # Add thresholds
    plt.axhline(-np.log10(0.05), color='gray', linestyle='--', lw=1)
    plt.axvline(2, color='gray', linestyle='--', lw=1)
    plt.axvline(-2, color='gray', linestyle='--', lw=1)

    plt.xlim(-6, 4.5)
    plt.ylim(0, 4.5)
    plt.title(title)
    plt.xlabel('logFC')
    plt.ylabel('-log10(adjusted p-value)')
    plt.tight_layout()
    plt.savefig(save_as)
    plt.close()

    # Gray-only version
    gray_save = save_as.replace('.svg', '_gray.svg')
    plt.figure(figsize=(6, 8))
    sns.scatterplot(data=df, x='logFC', y='-log10(padj)',
                    color='gray', alpha=0.9, s=40, edgecolor=None, legend=False)
    plt.axhline(-np.log10(0.05), color='gray', linestyle='--', lw=1)
    plt.axvline(1, color='gray', linestyle='--', lw=1)
    plt.axvline(-1, color='gray', linestyle='--', lw=1)
    plt.xlim(-6, 4.5)
    plt.ylim(0, 4.5)
    plt.title(title + " (gray only)")
    plt.xlabel('logFC')
    plt.ylabel('-log10(adjusted p-value)')
    plt.tight_layout()
    plt.savefig(gray_save)
    plt.close()

    # Return top 5 colored significant taxa
    top = [t for t, _ in sorted(
        ((tax, r['adj.P.Val']) for tax, r in df.iterrows() if r.color!='gray' and r['-log10(padj)']>2),
        key=lambda x: x[1]
    )[:5]]
    return top

def plot_correlations(inoculum_da_taxa,
                      drought_signature,
                      ds,
                      save_as):
    all_correlations = []
    for inoculum_da in inoculum_da_taxa:
        for signature_taxa in drought_signature:
            if (signature_taxa in ds.get_counts_features() and inoculum_da not in signature_taxa):
                res = spearmanr(ds.get_counts_feature(signature_taxa),
                                ds.get_counts_feature(inoculum_da))
                if res.pvalue <= 0.05:
                    all_correlations.append(res.statistic)

    if all_correlations:
        plt.figure(figsize=(8, 6))
        sns.kdeplot(all_correlations, color="dimgray", fill=True, alpha=0.7)
        plt.xlabel("Spearman")
        plt.ylabel("Density")
        plt.title(f"Distribution of significant Spearman correlations\n({len(all_correlations)}")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(save_as)
        plt.close()
    else:
        print("[WARNING] No significant correlations found.")

# node_colors maps each node name to a color
def draw_network(edge_df,
                 node_colors,
                 save_as):
    G = nx.from_pandas_edgelist(edge_df,
                                source='source',
                                target='target',
                                edge_attr='weight')

    colors = [node_colors.get(node, 'gray') for node in G.nodes()]

    edge_weights = nx.get_edge_attributes(G, 'weight')
    edge_colors = ['green' if w > 0 else 'red' for w in edge_weights.values()]
    edge_widths = [abs(w) for w in edge_weights.values()]

    pos = nx.random_layout(G, seed=42)

    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=300)

    nx.draw_networkx_edges(G, pos,
                           edge_color=edge_colors,
                           width=edge_widths)

    nx.draw_networkx_labels(G, pos, font_size=10)

    plt.axis('off')
    plt.tight_layout()
    plt.savefig(save_as)
    plt.close()
