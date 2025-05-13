import yaml
import os
import argparse

from preprocess_and_filter import preprocess_and_filter, process_metadata
from merge_datasets import merge_datasets
from plotting import *
from analysis import process_permanova, run_limma_diff_abundance_treatment, run_external_signature_validation

from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import LabelEncoder

from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_predict
from scipy.stats import spearmanr, pearsonr
import numpy as np
from sklearn.pipeline import Pipeline

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import networkx as nx

from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.distance import permanova
from scipy.spatial.distance import pdist, squareform

with open('config.yml') as f:
    config = yaml.safe_load(f)

def merge_stats():
    for level, level_name in config["levels"].items():
        feature_stats_outfile = os.path.join(config["data_path"], f"feature_stats_l{level}.tsv")
        sample_stats_outfile = os.path.join(config["data_path"], f"sample_stats_l{level}.tsv")
        feature_lines = []
        sample_lines = []
        header_written_feature = False
        header_written_sample = False
        for study in config["drought_studies"]:
            study_path = os.path.join(config["data_path"], study)
            feature_file = os.path.join(study_path,
                                        f"feature_filtering_stats_l{level}.tsv")
            sample_file = os.path.join(study_path,
                                       f"sample_filtering_stats_l{level}.tsv")
            if os.path.exists(feature_file):
                with open(feature_file, 'r') as f:
                    lines = f.readlines()
                if not header_written_feature:
                    feature_lines.append(lines[0].strip())
                    header_written_feature = True
                feature_lines.append(lines[1].strip())
            if os.path.exists(sample_file):
                with open(sample_file, 'r') as f:
                    lines = f.readlines()
                if not header_written_sample:
                    sample_lines.append(lines[0].strip())
                    header_written_sample = True
                sample_lines.append(lines[1].strip())
        with open(feature_stats_outfile, 'w') as f:
            f.write("\n".join(feature_lines) + "\n")
        with open(sample_stats_outfile, 'w') as f:
            f.write("\n".join(sample_lines) + "\n")

def compute_pairwise_lfc(counts_df, meta_df, group_col, group1, group2, taxon, max_pairs=500):
    # Ensure taxon is a column
    if taxon not in counts_df.columns:
        return np.array([])

    # Compute CPM per sample, log2-transform
    cpm = counts_df.div(counts_df.sum(axis=1), axis=0) * 1e6
    logcpm = np.log2(cpm + 1)

    # Samples in each group
    g1 = meta_df.index[meta_df[group_col] == group1]
    g2 = meta_df.index[meta_df[group_col] == group2]
    g1 = [s for s in g1 if s in logcpm.index]
    g2 = [s for s in g2 if s in logcpm.index]
    if not g1 or not g2:
        return np.array([])

    # Sample pairs
    pairs = min(max_pairs, len(g1) * len(g2))
    idx1 = np.random.choice(g1, size=pairs, replace=True)
    idx2 = np.random.choice(g2, size=pairs, replace=True)

    # Compute LFC for taxon by pairing samples
    return logcpm.loc[idx1, taxon].values - logcpm.loc[idx2, taxon].values

def compute_alpha_beta_bars(df, label):
    # Separate metadata and count data
    meta = df[['Inoculum']].copy()
    counts = df.drop(columns=[col for col in df.columns if col not in meta.columns])

    # Keep only numeric (taxa) columns
    counts = counts.select_dtypes(include=[np.number])

    # Compute alpha diversity (Shannon)
    alpha = alpha_diversity('shannon', counts.values, ids=counts.index)

    # Compute beta diversity (Bray-Curtis)
    beta_dm = beta_diversity('braycurtis', counts.values, ids=counts.index)

    # Compute within-group average distances
    results = []
    for group in ['DroughtLegacy', 'Other']:
        sample_ids = meta[meta['Inoculum'] == group].index
        sample_ids = [s for s in sample_ids if s in counts.index]

        if len(sample_ids) < 2:
            continue

        # Alpha: average per sample
        alpha_vals = alpha.loc[sample_ids]
        alpha_mean = alpha_vals.mean()

        # Beta: average pairwise distance within group
        sub_dm = beta_dm.filter(sample_ids)
        dists = sub_dm.condensed_form()
        beta_mean = dists.mean()

        results.append({
            'Inoculum': group,
            'Dataset': label,
            'Alpha Diversity': alpha_mean,
            'Beta Diversity': beta_mean
        })

    return pd.DataFrame(results)

def run_preprocessing():
    processed_files = dict()

    processed_files[0] = dict()
    processed_files[1] = dict()
    for level in config["levels"]:
        print(f"[INFO] Level: {level}")
        processed_files[0][level] = dict()
        processed_files[1][level] = dict()
        
        for study in config["drought_studies"]:
            print(f"[INFO] Study: {study}")
            study_path = os.path.join(config["data_path"],
                                      study)
            processed_tsv, processed_tsv_counts = preprocess_and_filter(study_path,
                                                                        study,
                                                                        level,
                                                                        config["levels"][level][0] + "__")
            processed_files[0][level][study] = processed_tsv
            processed_files[1][level][study] = processed_tsv_counts
    
    with open(config["processed_files"], 'w') as f:
        yaml.dump(processed_files, f, default_flow_style=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preprocess", action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    preprocess = args.preprocess

    if (not os.path.exists(config["processed_files"])) or preprocess:
        run_preprocessing()

    with open(config["processed_files"]) as f:
        processed_files = yaml.safe_load(f)

    metadata_files = dict()
    studies = config["drought_studies"]

    for study in studies:
        print(study)
        study_path = os.path.join(config["data_path"], study)
        metadata_file = process_metadata(study_path,
                                         study,
                                         config[study])
        metadata_files[study] = metadata_file

    for level, level_name in config["levels"].items():
        merge_stats()

    
    for level, level_name in [[6, "genus"]]:#config["levels"].items():
        merged_file_normalized = merge_datasets(config["data_path"],
                                                level,
                                                processed_files[0][level],
                                                metadata_files,
                                                "normalized")
        merged_file_counts = merge_datasets(config["data_path"],
                                            level,
                                            processed_files[1][level],
                                            metadata_files,
                                            "counts")

        process_permanova(level,
                          config)

        # Differential abundance analysis
        diff_abundance_results = run_limma_diff_abundance_treatment(level,
                                                                    merged_file_counts)
        # Filter out "uncultured" bacteria
        # (this should have been done at the beginning so I need to fix this)
        filtered_results = {key: diff_abundance_results[key]
                            for key in diff_abundance_results
                            if "uncultured" not in key}
        diff_abundance_results = filtered_results

        diff_abundance_results_inoculum = dict()
        success = ["moore2023microbial", "zhang2022cross"]
        fail = ["swift2024drought", "munoz-ucros2021drought"]
        outlier_taxa = dict()
        all_outlier_taxa = dict()
        for studies, outcome in [[success, "_success"], [fail, "_fail"]]:
            processed_inoculum_files_counts = dict()
            metadata_inoculum_files = dict()
            
            for study in studies:
                study_path = os.path.join(config["data_path"],
                                          study)
                _, filtered_tsv_counts = preprocess_and_filter(study_path,
                                                               study,
                                                               level,
                                                               config["levels"][level][0] + "__")
                processed_inoculum_files_counts[study] = filtered_tsv_counts

                # process metadata
                metadata_file = process_metadata(study_path,
                                                 study,
                                                 config[study])
                metadata_inoculum_files[study] = metadata_file

            merged_file_counts_inoculum = merge_datasets(config["data_path"],
                                                         level,
                                                         processed_inoculum_files_counts,
                                                         metadata_inoculum_files,
                                                         "counts" + outcome,
                                                         inoculum_studies=True)
            model = "Inoculum + Study"
            diff_abundance_results_inoculum[outcome] = run_limma_diff_abundance_treatment(level,
                                                                                          merged_file_counts_inoculum,
                                                                                          p_val=1,
                                                                                          formula=model,
                                                                                          variable="Inoculum",
                                                                                          base="Other",
                                                                                          condition="DroughtLegacy")
            save_as = os.path.join(config["plotting_dir"],
                                   f"volcano_plot{outcome}.svg")
            all_outlier_taxa[outcome] = make_volcano_plot(diff_abundance_results_inoculum[outcome],
                                                          diff_abundance_results,
                                                          save_as)
            outlier_taxa[outcome] = all_outlier_taxa[outcome][:5]

        print(outlier_taxa)
        plot_lfc_diff_abundance(diff_abundance_results,
                                level,
                                config)
        


    # Load dataframes that already contain metadata in columns
    original_df = pd.read_csv(
        os.path.join(config['data_path'], 'merged_taxonomy_counts_l6.tsv'), sep='\t', index_col=0
    )
    success_df = pd.read_csv(
        os.path.join(config['data_path'], 'merged_taxonomy_counts_success_l6.tsv'), sep='\t', index_col=0
    )
    fail_df = pd.read_csv(
        os.path.join(config['data_path'], 'merged_taxonomy_counts_fail_l6.tsv'), sep='\t', index_col=0
    )

    dfs = {'Success': success_df, 'Fail': fail_df, 'Original': original_df}
    grouping = {
        'Original': ('Treatment', 'Drought', 'Control'),
        'Success': ('Inoculum', 'DroughtLegacy', 'Other'),
        'Fail': ('Inoculum', 'DroughtLegacy', 'Other')
    }

    # Combine all outlier taxa
    taxa_all = list({
        *outlier_taxa.get('_success', []),
        *outlier_taxa.get('_fail', []),
        *outlier_taxa.get('original', [])
    })

    # Prepare long-form DataFrame for violin plot
    plot_rows = []
    for label, df_full in dfs.items():
        group_col, g1, g2 = grouping[label]
        meta_df = df_full[[group_col]]
        counts_df = df_full.drop(columns=[col for col in df_full if "p__" not in col])
        for taxon in taxa_all:
            lfc_vals = compute_pairwise_lfc(counts_df, meta_df, group_col, g1, g2, taxon)
            for val in lfc_vals:
                plot_rows.append({'Taxon': taxon, 'Dataset': label, 'LFC': val})

    plot_df = pd.DataFrame(plot_rows)

    # Plot violin: for each Taxon, Datasets side-by-side
    plt.figure(figsize=(max(6, len(taxa_all) * 1.5), 6))
    sns.violinplot(
        data=plot_df,
        x='Taxon',
        y='LFC',
        hue='Dataset',
        order=taxa_all,
        hue_order=['Success', 'Fail', 'Original'],
        dodge=True,
        cut=0
    )
    plt.axhline(0, color='gray', linestyle='--')
    plt.xlabel('Taxon')
    plt.ylabel('log2 fold change')
    plt.title('Pairwise log2-Fold Change Distributions')
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Dataset', loc='upper right')
    plt.tight_layout()

    out_png = os.path.join(config['plotting_dir'], 'violin_lfc_all_taxa.png')
    plt.savefig(out_png)
    plt.close()

    # Run for both datasets
    alpha_beta_df = pd.concat([
        compute_alpha_beta_bars(success_df, 'Success'),
        compute_alpha_beta_bars(fail_df, 'Fail')
    ])

    # Convert to long-form for plotting
    alpha_long = alpha_beta_df.melt(
        id_vars=['Inoculum', 'Dataset'],
        value_vars=['Alpha Diversity', 'Beta Diversity'],
        var_name='Metric',
        value_name='Diversity'
    )

    # Plot
    plt.figure(figsize=(8, 6))
    sns.barplot(
        data=alpha_long,
        x='Inoculum',
        y='Diversity',
        hue='Metric',
        palette='muted',
        ci='sd',
        dodge=True
    )
    plt.title("Alpha and Beta Diversity by Inoculum Type")
    plt.xlabel("Inoculum")
    plt.ylabel("Diversity Index")
    plt.legend(title='Diversity Metric')
    plt.tight_layout()

    plt.savefig(os.path.join(config["plotting_dir"], 'alpha_beta_diversity_by_inoculum.png'))
    plt.close()


if __name__ == "__main__":
    main()
