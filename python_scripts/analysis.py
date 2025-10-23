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
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def run_maaslin_diff_abundance(dataset,
                               fixed_effects="Treatment",
                               random_effects=None,
                               normalization="TSS",
                               transform="LOG",
                               correction="BH",
                               min_abundance=0.0,
                               min_prevalence=0.0,
                               p_val=1.0,
                               base="Control",
                               condition="Drought"):
    pandas2ri.activate()
    maaslin2 = importr('Maaslin2')

    counts_df = dataset.get_counts_features()
    metadata_df = dataset.get_metadata()

    counts_r = pandas2ri.py2rpy(counts_df.T)
    metadata_r = pandas2ri.py2rpy(metadata_df)

    robjects.globalenv['v'] = counts_r
    robjects.globalenv['metadata'] = metadata_r

    random_effects_str = f'c("{random_effects}")' if random_effects else 'NULL'

    robjects.r(f'''
    library(Maaslin2)

    metadata${fixed_effects} <- as.factor(metadata${fixed_effects})
    metadata${fixed_effects} <- relevel(metadata${fixed_effects}, ref="{base}")

    fit_data <- Maaslin2(
        input_data = v,
        input_metadata = metadata,
        output = "maaslin2_out",
        fixed_effects = c("{fixed_effects}"),
        random_effects = {random_effects_str},
        normalization = "{normalization}",
        transform = "{transform}",
        min_abundance = {min_abundance},
        min_prevalence = {min_prevalence},
        correction = "{correction}",
        plot_scatter=FALSE
    )
    ''')

    results_r = robjects.r('fit_data$results')
    results_df = results_r

    condition_label = f"{fixed_effects}{condition}"
    results_df = results_df[
        (results_df['metadata'] == fixed_effects) &
        (results_df['value'] == condition)
    ]
    print(results_df)
    results_df = results_df[results_df['qval'] <= p_val]

    results = {}
    for _, row in results_df.iterrows():
        feature = row['feature']
        results[feature] = {
            'logFC': row['coef'],
            'adj.P.Val': row['qval']
        }

    print(f"[INFO] Found {len(results)} DA features with p={p_val} ({fixed_effects}: {condition} vs {base})")

    return results


# Based on the MMuPHin user tutorial: https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html
def run_mmuphin_diff_abundance(dataset):
    pandas2ri.activate()
        
    mmuphin = importr('MMUPHin')

    df = dataset.get_counts_features()
    metadata_df = dataset.get_metadata()
    # Samples are columns
    df_r = pandas2ri.py2rpy(df.T)
    # Samples are rows
    metadata_r = pandas2ri.py2rpy(metadata_df)
    robjects.globalenv['data'] = df_r
    robjects.globalenv['metadata'] = metadata_r

    robjects.r(f'''
    library(MMUPHin)
    library(magrittr)
    library(dplyr)
    fit_lm_meta <- lm_meta(feature_abd = data,
                           batch = "StudyID",
                           exposure = "Treatment",
                           data = metadata,
                           control = list(verbose = TRUE, transform="LOG"))
    
    meta_fits <- fit_lm_meta$meta_fits
    meta_fits_summary <- meta_fits %>% 
         select(feature, coef, qval.fdr)
    ''')

    meta_fits_summary_r = robjects.r['meta_fits_summary']
    meta_fits_summary_py = pandas2ri.rpy2py(meta_fits_summary_r)

    # Reformat
    results = {}
    for _, row in meta_fits_summary_py.iterrows():
        results[row['feature']] = {
            'logFC': row['coef'],          
            'adj.P.Val': row['qval.fdr']  
        }

    print(f"[INFO] Found {len(results)} DA features for {dataset.dataset_name} (MMUPHin)")
    return results

def run_wilcoxon_diff_abundance(dataset,
                                p_val=1,
                                variable="Treatment",
                                base="Control",
                                condition="Drought",
                                correction="fdr_bh"):
    counts_df = dataset.get_counts_features()
    metadata_df = dataset.get_metadata()

    base_samples = metadata_df.index[metadata_df[variable] == base]
    cond_samples = metadata_df.index[metadata_df[variable] == condition]

    results = []
    for feature in counts_df.columns:
        base_vals = counts_df.loc[base_samples, feature].values
        cond_vals = counts_df.loc[cond_samples, feature].values

        try:
            stat, pval = mannwhitneyu(base_vals, cond_vals, alternative="two-sided")
        except ValueError:
            pval = 1.0

        mean_base = np.mean(base_vals + 1)
        mean_cond = np.mean(cond_vals + 1)
        logfc = np.log2(mean_cond / mean_base)

        results.append({
            "feature": feature,
            "logFC": logfc,
            "pval": pval
        })

    results_df = pd.DataFrame(results)

    results_df["pval_adj"] = multipletests(results_df["pval"], method=correction)[1]
    results_df = results_df[results_df["pval_adj"] <= p_val]

    out = {row["feature"]: {
                "logFC": row["logFC"],
                "pval": row["pval"],
                "adj.P.Val": row["pval_adj"]
            }
           for _, row in results_df.iterrows()}

    print(f"[INFO] Found {len(out)} DA features for {dataset.dataset_name} with Wilcoxon p<={p_val}")
    return out

def process_permanova(level,
                      compartment,
                      config,
                      batch_corrected=False):
    suffix = "_batch_corrected" if batch_corrected else ""
    be_vars = ["StudyID"]
    bio_vars = ["HostSpecific", "Treatment"]

    r_squared = {}

    for var in be_vars + bio_vars:
        if compartment == "bulk_soil" and var == "HostSpecific": continue
        permanova_dir = os.path.join(
            config["data_path"],
            "permanova",
            f"{var}_l{level}_{compartment}{suffix}"
        )
        subdir = get_qiime_extract_dir(permanova_dir)
        tsv_file = os.path.join(subdir, "adonis.tsv")

        res = pd.read_csv(tsv_file, sep='\t', index_col=0)
        r2_val = res.loc[var, 'R2']
        p_val = res.loc[var, 'Pr(>F)']
        r_squared[var] = r2_val

        if p_val > 0.001:
            print(f"p-value for {var} is {p_val:.4f} > 0.001")

    plotting_dir = config["plotting_dir"]
    indiv_df = pd.DataFrame({'Variable': list(r_squared.keys()),
                             'R2': list(r_squared.values())})

    print(compartment, suffix)
    print(indiv_df["Variable"])
    print(indiv_df["R2"])
    plt.figure(figsize=(8, 4))
    sns.barplot(data=indiv_df,
                x='Variable',
                y='R2',
                color='black',
                dodge=False)
    plt.ylim(0, 0.5)
    plt.xticks(rotation=45, ha='right')
    plt.title(f'PERMANOVA R2 by Variable (Level {level}{suffix})')
    plt.tight_layout()
    plt.savefig(os.path.join(
        plotting_dir,
        f"permanova_plot_l{level}{suffix}_compartment.svg"
    ))
    plt.close()
    

