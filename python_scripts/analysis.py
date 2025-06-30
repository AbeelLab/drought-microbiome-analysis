import pandas as pd
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova
import matplotlib.pyplot as plt
import os
# trim_taxonomy should probably be in utils rather than plotting
from plotting import trim_taxonomy
import numpy as np
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from utils import get_qiime_extract_dir
import matplotlib.pyplot as plt
import seaborn as sns
from preprocess_and_filter import *
from itertools import combinations
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score

def process_permanova(level,
                      config):
    be_vars = ["Primers", "Study"]
    bio_vars = ["Host", "HostSpecific", "Treatment", "Location"]

    r_squared = {}

    for var in be_vars + bio_vars:
        permanova_dir = os.path.join(config["data_path"],
                                     "permanova",
                                     f"{var}_l{level}")
        subdir = get_qiime_extract_dir(permanova_dir)
        tsv_file = os.path.join(subdir, "adonis.tsv")

        res = pd.read_csv(tsv_file, sep='\t', index_col=0)
        r2_val = res.loc[var, 'R2']
        p_val = res.loc[var, 'Pr(>F)']
        r_squared[var] = r2_val

        # manual inspection indicates they are all 0.001 but just to be sure
        if p_val > 0.001:
            print(f"p-value for {var} is {p_val:.4f} > 0.001")

    plotting_dir = config["plotting_dir"]

    indiv_df = pd.DataFrame({'Variable': list(r_squared.keys()),
                             'R2': list(r_squared.values())})

    plt.figure(figsize=(8, 4))
    sns.barplot(data=indiv_df,
                x='Variable',
                y='R2',
                color='black',
                dodge=False)
    plt.ylim(0, 0.5)
    plt.xticks(rotation=45, ha='right')
    plt.title(f'PERMANOVA R2 by Variable (Level {level})')
    plt.tight_layout()
    plt.savefig(os.path.join(plotting_dir, f"permanova_plot_l{level}.svg"))
    plt.close()

    combinations = [
        "HostSpecific+Study",
        "Host+Study",
        "Host+HostSpecific+Treatment+Location+Study",
        "HostSpecific+Treatment+Location+Study",
        "Host+Treatment+Location+Study"
    ]

    for combo in combinations:
        vars_in_combo = combo.split('+')
        permanova_dir = os.path.join(config["data_path"],
                                     "permanova",
                                     f"{combo}_l{level}")
        subdir = get_qiime_extract_dir(permanova_dir)
        tsv_file = os.path.join(subdir, "adonis.tsv")

        res = pd.read_csv(tsv_file, sep='\t', index_col=0)
        values = {v: res.loc[v, 'R2'] for v in vars_in_combo}
        values['Residuals'] = res.loc['Residuals', 'R2']

        fig, ax = plt.subplots(figsize=(6, 1.5))
        left = 0
        for var, r2 in values.items():
            ax.barh(0, r2, left=left, label=var)
            left += r2
        ax.set_xlim(0, 1)
        ax.set_yticks([])
        ax.set_xlabel('R2')
        ax.set_title(f'Stacked PERMANOVA R2 (Level {level} - {combo})')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig(os.path.join(
            plotting_dir,
            f"permanova_stacked_l{level}_{combo.replace('+', '_')}.svg"
        ))
        plt.close(fig)

    return r_squared

def run_limma_diff_abundance(dataset,
                             p_val=0.05,
                             formula="Treatment + Study + HostSpecific + Study * HostSpecific",
                             variable="Treatment",
                             base="Control",
                             condition="Drought",
                             correction="BH"):
    pandas2ri.activate()

    # Load required R packages
    limma = importr('limma')
    edgeR = importr('edgeR')

    # limma input requires separate metadata and transposed abundance matrix
    counts_df = dataset.get_counts_features()
    metadata_df = dataset.get_metadata()

    assert not counts_df.isna().values.any(), "Count table includes NaNs"
    assert not metadata_df[variable].isna().values.any(), f"NaN {variable} values"
    
    counts_r = pandas2ri.py2rpy(counts_df.T)
    metadata_r = pandas2ri.py2rpy(metadata_df)
    robjects.globalenv['v'] = counts_r
    robjects.globalenv['metadata'] = metadata_r

    # Fit limma
    robjects.r(f'''
    library(edgeR)
    library(limma)

    dge <- DGEList(counts = v)
    dge <- calcNormFactors(dge)

    # Control is the baseline
    metadata${variable} <- as.factor(metadata${variable})
    metadata${variable} <- relevel(metadata${variable}, ref="{base}")
    design <- model.matrix(~ {formula}, data=metadata)

    # voom transformation
    v_voom <- voom(dge, design, plot=FALSE)

    # fit model
    fit <- lmFit(v_voom, design)
    fit <- eBayes(fit)

    tt <- topTable(fit, coef="{variable}{condition}", p.value={p_val}, number=Inf, adjust.method="{correction}")
    ''')

    tt_df = pandas2ri.rpy2py(robjects.globalenv['tt'])

    results = {}
    for feature, row in tt_df.iterrows():
        results[feature] = {
            'logFC': row['logFC'],
            'adj.P.Val': row['adj.P.Val']
        }

    print(f"[INFO] Found {len(results)} DA features for {dataset.dataset_name} with p={p_val}")

    return results
