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

def plot_top10_features_per_component(pca, trimmed_features, plotting_dir, level):
    comp = pca.components_
    n_components = comp.shape[0]
    for i in range(n_components):
        loadings = comp[i]
        abs_loadings = np.abs(loadings)
        top10_idx = np.argsort(abs_loadings)[::-1][:10]
        top10_features = trimmed_features[top10_idx]
        top10_values = abs_loadings[top10_idx]
        plt.figure(figsize=(8, 6))
        sns.barplot(x=top10_values, y=top10_features, orient='h', color="skyblue")
        plt.xlabel("Absolute feature importance")
        plt.ylabel("Feature")
        plt.title(f"Top 10 Features for PC{i+1} (Level {level})")
        plt.tight_layout()
        save_to = os.path.join(plotting_dir, f"pca_top10_features_PC{i+1}_l{level}.png")
        plt.savefig(save_to, dpi=600, bbox_inches="tight")
        print(f"[INFO] Top 10 features bar plot for PC{i+1} saved to {save_to}")
        plt.close()

def plot_pca(level,
             level_name,
             merged_file,
             plotting_dir):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    taxonomic_features = df.columns.drop("Study")
    X = df[taxonomic_features]
    study_names = df["Study"]
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(X)
    var_exp = pca.explained_variance_ratio_
    comp = pca.components_
    trimmed_features = np.array([trim_taxonomy(f) for f in taxonomic_features])
    top_feature_pc1 = trimmed_features[np.argmax(np.abs(comp[0]))]
    top_feature_pc2 = trimmed_features[np.argmax(np.abs(comp[1]))]
    top_feature_pc3 = trimmed_features[np.argmax(np.abs(comp[2]))]
    top5_pc1_idx = np.argsort(np.abs(comp[0]))[::-1][:5]
    top5_pc2_idx = np.argsort(np.abs(comp[1]))[::-1][:5]
    top5_pc3_idx = np.argsort(np.abs(comp[2]))[::-1][:5]
    top5_pc1 = list(zip(trimmed_features[top5_pc1_idx], comp[0][top5_pc1_idx]))
    top5_pc2 = list(zip(trimmed_features[top5_pc2_idx], comp[1][top5_pc2_idx]))
    top5_pc3 = list(zip(trimmed_features[top5_pc3_idx], comp[2][top5_pc3_idx]))
    print("Top 5 features for PC1:", top5_pc1)
    print("Top 5 features for PC2:", top5_pc2)
    print("Top 5 features for PC3:", top5_pc3)
    pca_df = pd.DataFrame(pca_result, columns=["PC1", "PC2", "PC3"], index=df.index)
    pca_df["Study"] = study_names
    def format_text(pc1, pc2, pc3, top1, top2):
        text = f"Variance explained: $\\mathbf{{PC1}}$: {pc1:.2%}, " \
               f"$\\mathbf{{PC2}}$: {pc2:.2%}, $\\mathbf{{PC3}}$: {pc3:.2%}.\n" \
               f"Top features: $\\mathbf{{PC1}}$: {top1}, $\\mathbf{{PC2}}$: {top2}"
        return "\n".join(textwrap.wrap(text, width=80))
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Study", palette="Set2", alpha=0.75)
    plt.title(f"PCA: PC1 vs PC2\nLevel: {level_name}")
    plt.legend(title="Study", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.subplots_adjust(right=0.8, bottom=0.35)
    plt.figtext(0.5, 0.10, format_text(var_exp[0], var_exp[1], var_exp[2], top_feature_pc1, top_feature_pc2),
                ha="center", fontsize=10, wrap=True)
    save_to = os.path.join(plotting_dir, f"pca_plot_PC1_PC2_l{level}.png")
    plt.savefig(save_to, dpi=600, bbox_inches="tight")
    print(f"[INFO] PCA plot saved to {save_to}")
    plt.close()
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=pca_df, x="PC2", y="PC3", hue="Study", palette="Set2", alpha=0.75)
    plt.title(f"PCA: PC2 vs PC3\nLevel: {level_name}")
    plt.legend(title="Study", bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.subplots_adjust(right=0.8, bottom=0.35)
    plt.figtext(0.5, 0.10, format_text(var_exp[0], var_exp[1], var_exp[2], top_feature_pc2, top_feature_pc3),
                ha="center", fontsize=10, wrap=True)
    save_to = os.path.join(plotting_dir, f"pca_plot_PC2_PC3_l{level}.png")
    plt.savefig(save_to, dpi=600, bbox_inches="tight")
    print(f"[INFO] PCA plot saved to {save_to}")
    plt.close()
    plot_top10_features_per_component(pca, trimmed_features, plotting_dir, level)

def sparsity_heatmap(level,
                     level_name,
                     merged_file,
                     plotting_dir,
                     log_scale=False,
                     shuffle=False):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    taxonomic_features_df = df.drop(columns=["Study"])

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
