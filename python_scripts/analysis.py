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

def run_limma_diff_abundance_treatment(level,
                                       merged_file,
                                       p_val=0.001,
                                       formula="Treatment + Study + HostSpecific",
                                       variable="Treatment",
                                       base="Control",
                                       condition="Drought"):
    df = pd.read_csv(merged_file, sep='\t', index_col=0)
    print(df)

    pandas2ri.activate()

    # Load required R packages
    limma = importr('limma')
    edgeR = importr('edgeR')

    # Dependent variables
    taxonomic_features = [col for col in df.columns if "p__" in col]

    # Variables (columns) that may lead to undesired batch effects
    # Because Primers and regions are so dependent on Study, they should be removed
    # if batch effects for Study are removed
    be_vars = ["Study"]

    # Variables that lead to biological variation (should be kept)
    bio_vars = ["HostSpecific", "Location", "SoilType", "Treatment", "Inoculum"]

    # limma input requires separate metadata and transposed abundance matrix
    abundance_df = df[taxonomic_features]
    metadata_df = df[bio_vars + be_vars].astype("category")
    metadata_df = metadata_df.dropna(subset=['Treatment'])
    abundance_df = abundance_df.loc[metadata_df.index, :]
    abundance_r = pandas2ri.py2rpy(abundance_df.T)
    metadata_r = pandas2ri.py2rpy(metadata_df)
    robjects.globalenv['v'] = abundance_r
    robjects.globalenv['metadata'] = metadata_r

    # Fit limma
    robjects.r(f'''
    library(edgeR)
    library(limma)

    dge <- DGEList(counts = v)
    dge <- calcNormFactors(dge)

    # Control is the baseline
    metadata${variable} <- relevel(metadata${variable}, ref="{base}")
    design <- model.matrix(~ {formula}, data=metadata)

    # voom transformation
    v_voom <- voom(dge, design, plot=FALSE)

    # fit model
    fit <- lmFit(v_voom, design)
    fit <- eBayes(fit)

    tt <- topTable(fit, coef="{variable}{condition}", p.value={p_val}, number=Inf, adjust.method="BH")
    ''')

    tt_df = pandas2ri.rpy2py(robjects.globalenv['tt'])

    results = {}
    for feature, row in tt_df.iterrows():
        results[feature] = {
            'logFC': row['logFC'],
            'adj.P.Val': row['adj.P.Val']
        }

    return results


def run_external_signature_validation(diff_abundance_results,
                                      config,
                                      level=6,
                                      level_shortcut="g__"):
    unfiltered_dfs = dict()
    filtered_dfs = dict()
    metadata_dfs = dict()
    for study in config["inoculum_studies"]:
        study_path = os.path.join(config["data_path"],
                                  study)
        # filter features based on abundance/prevalence
        # this will also generate a file with all features (unfiltered)
        filtered_tsv, _ = preprocess_and_filter(study_path,
                                                study,
                                                level,
                                                level_shortcut)
        # I shouldn't do this
        unfiltered_tsv = filtered_tsv.replace("filtered",
                                              "before-feature-filtering")

        filtered_dfs[study] = pd.read_csv(filtered_tsv, sep='\t', index_col=0).drop('Study', axis=1)
        unfiltered_dfs[study] = pd.read_csv(unfiltered_tsv, sep='\t', index_col=0).drop("Study", axis=1)

        # process metadata
        metadata_file = process_metadata(study_path,
                                         study,
                                         config[study])
        metadata_dfs[study] = pd.read_csv(metadata_file, sep='\t', index_col=0)

    # intersect features from all studies
    all_features_intersection = set.intersection(*[set(df.columns.tolist())
                                                   for df in unfiltered_dfs.values()])
    filtered_features_intersection = list(set.intersection(*[set(df.columns.tolist())
                                                             for df in filtered_dfs.values()]))

    # intersect features with signature
    signature_taxa = diff_abundance_results.keys()
    signature_intersection = list(all_features_intersection.intersection(signature_taxa))
    all_features_intersection = list(all_features_intersection)

    # reduce dataset sizes to intersection
    unfiltered_dfs = {study: df[all_features_intersection]
                      for study, df in unfiltered_dfs.items()}
    filtered_dfs = {study: df[filtered_features_intersection]
                    for study, df in filtered_dfs.items()}
    signature_dfs = {study: df[signature_intersection]
                    for study, df in unfiltered_dfs.items()}
    # normalize signature features
    signature_dfs = {study: df.div(df.sum(axis=1), axis=0) * 100
                    for study, df in signature_dfs.items()}

    experiments = {"Unfiltered feature set": unfiltered_dfs,
                   "Filtered feature set": filtered_dfs,
                   "Signature feature set": signature_dfs}

    print(unfiltered_dfs['swift2024drought'])
    print(filtered_dfs['swift2024drought'])
    print(signature_dfs['swift2024drought'])

    performance = dict()
    for exp_name, dfs in experiments.items():
        performance[exp_name] = {}
        studies = list(dfs.keys())

        for test_study in studies:
            train_studies = [s for s in studies if s != test_study]
            X_train = pd.concat([dfs[s] for s in train_studies], axis=0)
            # predict Treatment (Drought=1 vs Control=0)
            y_train = pd.concat([metadata_dfs[s].loc[dfs[s].index,
                                                     'Treatment'].apply(lambda x: 1 if x == 'Drought' else 0)
                                 for s in train_studies])
            # test study is the one left over
            X_test = dfs[test_study]
            y_test = metadata_dfs[test_study].loc[X_test.index,
                                                  'Treatment'].apply(lambda x: 1 if x == 'Drought' else 0)

            # fit RF
            clf = RandomForestClassifier(n_estimators=100, random_state=0)
            clf.fit(X_train.values, y_train.values)
            preds = clf.predict(X_test.values)
            performance[exp_name][test_study] = accuracy_score(y_test.values, preds)

    print(performance)
    performance_df = pd.DataFrame(performance).T
    print(performance_df)

    # x axis: test study name
    # y axis: feature set (number of features)
    
    return performance_df
