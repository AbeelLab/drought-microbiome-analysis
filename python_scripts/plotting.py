import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import LogNorm
import seaborn as sns
from sklearn.decomposition import PCA
import os
import numpy as np
import textwrap

def trim_taxonomy(tax_str):
    parts = tax_str.split(";")
    for part in reversed(parts):
        part = part.strip()
        if part not in ["", "__"]:
            return part
    return ""

def plot_pca(level,
             level_name,
             merged_file,
             plotting_dir):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)

    metadata = ["Study", "Treatment", "Host", "Host (specific)", "Inoculum"]
    taxonomic_features = df.columns.drop(metadata)
    X = df[taxonomic_features]
    other_features = df[metadata]

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(X)
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2"], index=df.index)
    pca_df = pd.concat([pca_df, other_features], axis=1)

    def save_pca_plot(color_by):
        plt.figure(figsize=(8, 6))
        
        sns.scatterplot(data=pca_df,
                        x="PC1", y="PC2",
                        hue=color_by,
                        palette="Paired",
                        alpha=0.75)
        plt.title(f"PCA: PC1 vs PC2\nLevel: {level_name} (Colored by {color_by})")
        
        plt.legend(title=color_by, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        save_to = os.path.join(plotting_dir, f"pca_plot_PC1_PC2_l{level}_by_{color_by.lower()}.png")
        plt.savefig(save_to, dpi=600, bbox_inches="tight")
        print(f"[INFO] PCA plot saved to {save_to}")
        plt.close()

    for category in ["Study", "Treatment", "Host", "Host (specific)"]:
        save_pca_plot(category)
    

def sparsity_heatmap(level,
                     level_name,
                     merged_file,
                     plotting_dir,
                     log_scale=False,
                     shuffle=False):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    taxonomic_features_df = df.drop(columns=["Study",
                                             "Treatment",
                                             "Host",
                                             "Inoculum",
                                             "Host (specific)"])

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
    
    # Create a grouping column: "Study | Host"
    df['Group'] = df['Study'] + " | " + df['Host']
    
    # Count samples per group and treatment
    count_df = df.groupby(['Group', 'Treatment']).size().unstack(fill_value=0)
    count_df = count_df.reset_index()
    count_df[['Study', 'Host']] = count_df['Group'].str.split(r' \| ', expand=True)
    count_df = count_df.sort_values(by=['Host', 'Study'])
    
    # Determine group order and positions
    groups = count_df['Group']
    y_positions = range(len(groups))
    
    # Get unique hosts and assign colors from the Paired palette
    unique_hosts = count_df['Host'].unique()
    host_palette = dict(zip(unique_hosts, sns.color_palette("Set2", len(unique_hosts))))
    
    # Get counts for treatments, ensuring both columns exist
    control_counts = count_df.get('Control', pd.Series([0] * len(count_df)))
    drought_counts = count_df.get('Drought', pd.Series([0] * len(count_df)))
    
    # Create list of base colors based on host for each group
    base_colors = [host_palette[host] for host in count_df['Host']]
    
    fig, ax = plt.subplots(figsize=(10, 0.5 * len(groups)))
    
    # Plot stacked bars: lower segment for Control (lower alpha), upper for Drought (full alpha)
    ax.barh(y_positions, control_counts, color=[(r, g, b, 0.5) for r, g, b in base_colors],
            edgecolor='black', label='Control')
    ax.barh(y_positions, drought_counts, left=control_counts,
            color=[(r, g, b, 1.0) for r, g, b in base_colors],
            edgecolor='black', label='Drought')
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(groups)
    ax.set_xlabel("Number of Samples")
    ax.set_title("Number of Samples per Study | Host")
    
    # Compute overall sample statistics
    total_samples = len(df)
    total_drought = (df['Treatment'] == 'Drought').sum()
    total_control = (df['Treatment'] == 'Control').sum()
    annotation = f"Total Samples: {total_samples} | Drought: {total_drought} | Control: {total_control}"
    
    plt.figtext(0.5, 0.01, annotation, wrap=True, horizontalalignment='center', fontsize=10)
    
    filename = "samples_stats.png"
    save_to = os.path.join(plotting_dir, filename)
    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(save_to, dpi=600)
    print(f"[INFO] Sample stats plot saved to {save_to}")
    plt.close()


def plot_taxonomy(level, merged_file, plotting_dir):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    df['Group'] = df['Study'] + " | " + df['Host']
    
    metadata = ["Study", "Treatment", "Host", "Host (specific)", "Inoculum"]
    taxonomic_features = df.columns.difference(metadata + ['Group'])
    
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
