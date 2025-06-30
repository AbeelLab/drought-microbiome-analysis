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

def trim_taxonomy(tax_str):
    parts = tax_str.split(";")
    for part in reversed(parts):
        part = part.strip()
        if part not in ["", "__"]:
            return part
    return ""

def plot_pcoa(level,
              level_name,
              merged_file,
              config):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)

    taxonomic_features = [col for col in df.columns
                          if "p__" in col]
    X = df[taxonomic_features]
    other_features = df[[col for col in df.columns
                         if "p__" not in col]]

    bc_distances = pdist(X.values, metric='braycurtis')
    distance_matrix = DistanceMatrix(squareform(bc_distances), ids=df.index)

    pcoa_object = pcoa(distance_matrix)
    pcoa_result = pcoa_object.samples[["PC1", "PC2"]]
    pcoa_df = pd.DataFrame(pcoa_result,
                           columns=["PC1", "PC2"],
                           index=df.index)
    pcoa_df = pd.concat([pcoa_df, other_features], axis=1)

    plotting_dir = config["plotting_dir"]
    def save_pcoa_plot(color_by):
        plt.figure(figsize=(8, 6))
        
        sns.scatterplot(data=pcoa_df,
                        x="PC1", y="PC2",
                        hue=color_by,
                        palette=config["colors"],
                        alpha=0.75)
        plt.title(f"PCoA: PC1 vs PC2\nLevel: {level_name} (legend: {color_by})\nBray-Curtis dissimilarity")
        
        plt.legend(title=color_by,
                   bbox_to_anchor=(1.05, 1),
                   loc='upper left',
                   borderaxespad=0)
        save_to = os.path.join(plotting_dir, f"pcoa_plot_PC1_PC2_l{level}_by_{color_by.lower()}.png")
        plt.savefig(save_to, dpi=600, bbox_inches="tight")
        print(f"[INFO] PCoA plot saved to {save_to}")
        plt.close()

    categories = ["Study (full name)", "Host"]

    for category in categories:
        save_pcoa_plot(category)
    

def sparsity_heatmap(level,
                     level_name,
                     merged_file,
                     plotting_dir,
                     log_scale=False,
                     shuffle=False):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    taxonomic_features = [col for col in df.columns
                          if "p__" in col]
    taxonomic_features_df = df[taxonomic_features]

    if shuffle:
        taxonomic_features_df = taxonomic_features_df.sample(frac=1, axis=1)
        shuffle_add_on = "_shuffled"
    else:
        shuffle_add_on = ""
        
    total_entries = taxonomic_features_df.size
    zero_entries = ((taxonomic_features_df == 0)
                    .astype(int)
                    .sum(axis=1)
                    .sum())
    sparsity = zero_entries / total_entries * 100

    cmap = sns.color_palette("Greys", as_cmap=True)
    if log_scale:
        log_add_on = "_log"
        smallest_nonzero = (taxonomic_features_df[taxonomic_features_df > 0]
                            .min()
                            .min())
        epsilon = smallest_nonzero / 10
        data = taxonomic_features_df.replace(0, epsilon)
        norm = LogNorm(vmin=epsilon, vmax=100)
    else:
        log_add_on = ""
        data = taxonomic_features_df
        norm = None

    ax = sns.heatmap(data,
                     cmap=cmap,
                     vmin=0 if not log_scale else epsilon,
                     vmax=100,
                     norm=norm,
                     cbar_kws={'label': 'Abundance (0 to 100%)'})

    rect = patches.Rectangle((0, 0), 1, 1,
                             transform=ax.transAxes,
                             fill=False,
                             color="black",
                             linewidth=2)
    ax.add_patch(rect)

    ax.set_title(f"Taxonomic features across samples and studies\n"
                 f"Level: {level_name}\n"
                 f"Sparsity: {sparsity:.2f}% zeros",
                 fontsize=10)
    ax.tick_params(left=False, bottom=False,
                   labelleft=False, labelbottom=False)
    ax.set_xlabel(f"Taxonomic features: {len(taxonomic_features_df.columns)}")
    ax.set_ylabel(f"Samples: {len(taxonomic_features_df)}")

    darkest_col = taxonomic_features_df.sum(axis=0).idxmax()
    wrapped_darkest = "\n".join(part.strip() for part in darkest_col.split(";"))
    ax.text(0.5, -0.12,
            f"Darkest column:\n{wrapped_darkest}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=10)

    plt.tight_layout()
    filename = f"sparsity_plot_l{level}{log_add_on}{shuffle_add_on}.png"
    save_to = os.path.join(plotting_dir, filename)
    plt.savefig(save_to, dpi=600)
    print(f"[INFO] Sparsity heatmap plot saved to {save_to}")
    plt.close()

def plot_num_samples(merged_file, plotting_dir):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    df['Group'] = df['Study (full name)'] + " | " + df['Host']
    count_df = df.groupby(['Group', 'Treatment']).size().unstack(fill_value=0)
    count_df = count_df.reset_index()
    count_df[['Study (full name)', 'Host']] = count_df['Group'].str.split(r' \| ', expand=True)
    count_df = count_df.sort_values(by=['Host', 'Study (full name)'])
    
    groups = count_df['Group']
    unique_hosts = count_df['Host'].unique()
    host_palette = dict(zip(unique_hosts, sns.color_palette("Set2", len(unique_hosts))))
    
    control_counts = count_df.get('Control', pd.Series([0] * len(count_df)))
    drought_counts = count_df.get('Drought', pd.Series([0] * len(count_df)))
    base_colors = [host_palette[host] for host in count_df['Host']]

    y_positions = []
    y_labels = []
    spacing = 0.5
    current_y = 0
    previous_host = None

    for i, (group, host) in enumerate(zip(groups, count_df['Host'])):
        if host != previous_host and previous_host is not None:
            current_y += spacing  # Extra gap between host groups
        y_positions.append(current_y)
        y_labels.append(group.split(" | ")[0])
        current_y += 1
        previous_host = host

    fig, ax = plt.subplots(figsize=(10, 0.5 * len(groups)))

    ax.barh(y_positions, control_counts, color=[(*rgb, 0.5) for rgb in base_colors],
            edgecolor='black', label='Control')
    ax.barh(y_positions, drought_counts, left=control_counts,
            color=[(*rgb, 1.0) for rgb in base_colors],
            edgecolor='black', label='Drought')

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels)
    ax.set_xlabel("Number of Samples")
    ax.set_title("Number of Samples per Study | Host")
    
    host_patches = [Patch(facecolor=host_palette[host], label=host) for host in unique_hosts]
    treatment_patches = [
        Patch(facecolor='grey', edgecolor='black', alpha=0.5, label='Control'),
        Patch(facecolor='grey', edgecolor='black', alpha=1.0, label='Drought')
    ]
    
    all_patches = treatment_patches + host_patches
    ax.legend(handles=all_patches, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    total_samples = len(df)
    total_drought = (df['Treatment'] == 'Drought').sum()
    total_control = (df['Treatment'] == 'Control').sum()
    annotation = f"Total Samples: {total_samples} | Drought: {total_drought} | Control: {total_control}"
    
    plt.figtext(0.5, 0.01, annotation, wrap=True, horizontalalignment='center', fontsize=10)
    
    filename = "samples_stats.png"
    save_to = os.path.join(plotting_dir, filename)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(save_to, dpi=600, bbox_inches='tight')
    plt.savefig(save_to.replace("png", "svg"))
    plt.close()

    # Pie chart: Host
    host_counts = df['Host'].value_counts()
    fig, ax = plt.subplots()
    ax.pie(host_counts, labels=host_counts.index, 
           colors=[host_palette[h] for h in host_counts.index],
           startangle=90, counterclock=False, autopct='%1.1f%%',
           wedgeprops={'edgecolor': 'white', 'linewidth': 1})
    ax.set_title("Samples per Host")
    save_host_pie = os.path.join(plotting_dir, "samples_per_host.png")
    plt.savefig(save_host_pie, dpi=600)
    plt.savefig(save_host_pie.replace("png", "svg"))
    plt.close()

    # Pie chart: Treatment
    treatment_counts = df['Treatment'].value_counts()
    treatment_colors = {'Control': '#d3d3d3', 'Drought': '#000000'}
    fig, ax = plt.subplots()
    ax.pie(treatment_counts, labels=treatment_counts.index,
           colors=[treatment_colors[t] for t in treatment_counts.index],
           startangle=90, counterclock=False, autopct='%1.1f%%',
           wedgeprops={'edgecolor': 'white', 'linewidth': 1})
    ax.set_title("Samples per Treatment")
    save_treatment_pie = os.path.join(plotting_dir, "samples_per_treatment.png")
    plt.savefig(save_treatment_pie, dpi=600)
    plt.savefig(save_treatment_pie.replace("png", "svg"))
    plt.close()

    print(f"[INFO] Sample stats plot saved to {save_to}")
    print(f"[INFO] Host pie chart saved to {save_host_pie}")
    print(f"[INFO] Treatment pie chart saved to {save_treatment_pie}")


def plot_taxonomy(level, merged_file, plotting_dir):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    df['Group'] = df['Study'] + " | " + df['Host']
    
    taxonomic_features = [col for col in df.columns
                          if "p__" in col]
    
    max_features = {2: 6, 3: 10, 4: 14, 5: 16, 6: 18}
    overall_means = df[taxonomic_features].mean(axis=0)
    top_features = overall_means.sort_values(ascending=False).head(max_features[level]).index.tolist()
    
    group_all = df.groupby('Group')[taxonomic_features].mean()
    total_all = group_all.sum(axis=1)
    
    group_top = df.groupby('Group')[top_features].mean()
    top_pct = group_top.div(total_all, axis=0) * 100
    top_pct = top_pct.fillna(0)
    
    other_pct = 100 - top_pct.sum(axis=1)
    other_pct[other_pct < 0] = 0
    
    group_pct = top_pct.copy()
    group_pct['Other'] = other_pct
    
    group_pct = group_pct.reset_index()
    group_pct[['Study', 'Host']] = group_pct['Group'].str.split(r' \| ', expand=True)
    group_pct = group_pct.sort_values(by=['Host', 'Study'])
    group_pct = group_pct.set_index('Group')
    
    groups = group_pct.index.tolist()
    y_positions = range(len(groups))
    
    fig, ax = plt.subplots(figsize=(10, 0.5 * len(groups)))
    
    feature_list = top_features + ['Other']
    palette = sns.color_palette("Paired", len(top_features))
    colors = dict(zip(top_features, palette))
    colors['Other'] = (0.8, 0.8, 0.8)
    
    cumulative = [0] * len(groups)
    for feature in feature_list:
        values = group_pct[feature].values
        ax.barh(y_positions, values, left=cumulative, color=colors[feature],
                edgecolor='black', label=feature)
        cumulative = [cum + val for cum, val in zip(cumulative, values)]

    ax.set_yticks(y_positions)
    ax.set_yticklabels(groups)
    ax.set_xlabel("Percentage (%)")
    ax.set_title("Taxonomic Composition (Top Features + Other)")
    ax.set_xlim(0, 100)
    
    # Prepare legend labels with last taxonomy rank
    def simplify_label(f):
        if f == 'Other':
            return 'Other'
        parts = f.split(';')
        non_empty = [p.split('__')[-1] for p in parts if p and p.split('__')[-1]]
        label = non_empty[-1] if non_empty else 'Unclassified'
        return "\n".join(label.split('-'))

    simplified_labels = [simplify_label(f) for f in feature_list]
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[f]) for f in feature_list]
    ax.legend(handles, simplified_labels, bbox_to_anchor=(1.05, 1), loc='upper left')
    
    filename = f"taxonomy_l{level}.png"
    save_to = os.path.join(plotting_dir, filename)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(save_to, dpi=600)
    print(f"[INFO] Sample stats plot saved to {save_to}")
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
