import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova
import matplotlib.pyplot as plt
import os

def run_permanova(level,
                  level_name,
                  merged_file,
                  plotting_dir=None):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    taxonomic_features = [col for col in df.columns if "p__" in col]
    metadata_columns_of_interest = [
        "Study", "Treatment", "Primers", "Host", "Location"
    ]

    tax_df = df[taxonomic_features]
    meta_df = df[metadata_columns_of_interest]

    dist_array = pdist(tax_df.values, metric="braycurtis")
    dist_matrix = squareform(dist_array)
    dm = DistanceMatrix(dist_matrix, ids=tax_df.index)

    stats = {}
    pvals = {}

    for col in metadata_columns_of_interest:
        col_df = meta_df[[col]].dropna()
        sub_dm = dm.filter(col_df.index)
        res = permanova(sub_dm, col_df, column=col, permutations=999)
        stats[col] = res["test statistic"]
        pvals[col] = res["p-value"]
        print(f"{col}: pseudo‑F = {res['test statistic']:.3f}, p‑value = {res['p-value']:.3f}")

    if plotting_dir:
        plt.figure(figsize=(8, 4))
        bars = plt.bar(stats.keys(), stats.values(), color='black', edgecolor='black')
        plt.xticks(rotation=45, ha='right')
        plt.ylabel("PERMANOVA pseudo‑F")
        plt.title(f"PERMANOVA pseudo‑F across metadata ({level_name})")

        for i, col in enumerate(stats.keys()):
            if round(pvals[col], 3) != 0.001:
                plt.text(i, stats[col] + 0.05, f"p={pvals[col]:.3f}",
                         ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        out_png = os.path.join(plotting_dir, f"permanova_{level_name}_pseudoF.png")
        plt.savefig(out_png, dpi=300)
        plt.close()
        print(f"[INFO] pseudo‑F barplot saved to {out_png}")

    return stats

def apply_batch_correction(level,
                           level_name,
                           original_file):
    df = pd.read_csv(original_file, sep='\t', index_col=0)

    # Dependent variables
    taxonomic_features = [col for col in df.columns if "p__" in col]

    # Variables (columns) that may lead to undesired batch effects
    # Because Primers and regions are so dependent on Study, they should be removed
    # if batch effects for Study are removed
    be_vars = ["Study"]

    # Variables that lead to biological variation (should be kept)
    bio_vars = ["Treatment", "Host (specific)", "Location", "Soil type"]

    # Apply batch correction with rpy limma, with bio_vars as covariates
    ...
    
    base, ext = os.path.splitext(original_file)
    corrected_file = f"{base}_batch_corrected{ext}"

    # Write df to corrected_file
    ...

    return corrected_file
