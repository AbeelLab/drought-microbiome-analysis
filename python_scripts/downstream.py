import yaml
import os
import argparse

from preprocess_and_filter import preprocess_and_filter, process_metadata
from merge_datasets import merge_datasets
from plotting import *
from analysis import process_permanova, run_limma_diff_abundance

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

        # print(merged_file)
        # Differential abundance analysis
        diff_abundance_results = run_limma_diff_abundance(level,
                                                          merged_file_counts)
        # Filter out "uncultured" bacteria
        # (this should have been done at the beginning so I need to fix this)
        filtered_results = {key: diff_abundance_results[key]
                        for key in diff_abundance_results
                        if "uncultured" not in key}
        diff_abundance_results = filtered_results
        plot_lfc_diff_abundance(diff_abundance_results,
                                level,
                                config)

        # TO DO: this should be cleaned up and moved to another file
        inoculum_diff_abundance_per_study = {}
        for study_name in config["inoculum_studies"]:
            # Process and filter
            study_path = os.path.join(config["data_path"],
                                      study_name)
            preprocess_and_filter(study_path,
                                  study_name,
                                  level,
                                  config["levels"][level][0] + "__")
            
            # Process metadata
            metadata_file = process_metadata(study_path,
                                             study_name,
                                             config[study_name])

            comp_dir = os.path.join(study_path, f"composition_table_l{level}")
            # Create 3 dataframes
            # One with unfiltered features (columns)
            f = os.path.join(comp_dir, "feature-table-before-feature-filtering.tsv")
            unfiltered_df = pd.read_csv(f, sep='\t', index_col=0)
            print("All features: ", len(unfiltered_df.columns))

            # One with filtered features (columns)
            f = os.path.join(comp_dir, "feature-table-filtered.tsv")
            filtered_df = pd.read_csv(f, sep='\t', index_col=0)
            print("Filtered features: ", len(filtered_df.columns))

            # One only with the features from among the leys of diff_abundance_results
            # (take feature intersection)
            feature_intersection = list(set(unfiltered_df.columns) & set(diff_abundance_results.keys()))
            print("Common: ", len(feature_intersection), " out of ", len(diff_abundance_results))
            signature_df = unfiltered_df[feature_intersection]
            # Re-normalize (TSS)
            signature_df = signature_df.div(signature_df.sum(axis=1), axis=0)

            metadata_df = pd.read_csv(metadata_file, sep='\t', index_col=0)
            unfiltered_df = unfiltered_df.join(metadata_df)
            filtered_df = filtered_df.join(metadata_df)
            signature_df = signature_df.join(metadata_df)

            unfiltered_df = unfiltered_df[unfiltered_df['Sample type'] == 'Plant-associated']
            filtered_df = filtered_df[filtered_df['Sample type']   == 'Plant-associated']
            signature_df = signature_df[signature_df['Sample type'] == 'Plant-associated']
            def eval_rf(X, y, name):
                cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=101)
                clf = RandomForestClassifier(random_state=101)
                scores = cross_val_score(
                    clf, X, y,
                    cv=cv,
                    scoring='accuracy',
                    n_jobs=-1
                )
                print(f"[{name}] accuracy: {scores.mean():.3f} ± {scores.std():.3f}")
                return scores

            # Evaluate on each DataFrame
            y = unfiltered_df["Treatment"]
            binarize = {"Control": 0, "Drought": 1}
            y = [binarize[label] for label in y]

            num_samples = len(y)
            print("# samples: ", num_samples)
            print("Drought samples: ", float(sum(y))/len(y))

            scores_unfiltered = eval_rf(unfiltered_df[[col for col in unfiltered_df.columns
                                                       if "p__" in col]],
                                        y,
                                        'Unfiltered')
            scores_filtered   = eval_rf(filtered_df[[col for col in filtered_df.columns
                                                     if "p__" in col]],
                                        y,
                                        'Filtered')
            scores_signature  = eval_rf(signature_df[[col for col in signature_df.columns
                                                      if "p__" in col]],
                                        y,
                                        'Signature')

            # Dummy classifier
            dummy = DummyClassifier(strategy='most_frequent')
            cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
            dummy_scores = cross_val_score(dummy,
                                           unfiltered_df[[col for col in unfiltered_df.columns
                                                          if "p__" in col]],
                                           y, cv=cv, scoring='accuracy', n_jobs=-1)
            print(f"Dummy accuracyc: {dummy_scores.mean():.3f} ± {dummy_scores.std():.3f}")

            # Prepare data for plotting
            data = {
                "Accuracy": list(scores_unfiltered) + list(scores_filtered) + list(scores_signature),
                "Classifier": (["All"] * len(scores_unfiltered) +
                               ["Filtered"] * len(scores_filtered) +
                               ["Drought signature"] * len(scores_signature))}
            plot_df = pd.DataFrame(data)

            n_unfiltered = sum("p__" in col for col in unfiltered_df.columns)
            n_filtered   = sum("p__" in col for col in filtered_df.columns)
            n_signature  = sum("p__" in col for col in signature_df.columns)

            plt.figure(figsize=(3, 6))
            sns.boxplot(data=plot_df, x="Classifier", y="Accuracy", color="#CCCCCC")

            # Dummy score as horizontal line
            dummy_median = np.median(dummy_scores)
            plt.axhline(dummy_median, color="black", linestyle="--", label=f"Dummy median = {dummy_median:.2f}")

            xtick_labels = [
                f"All\n({n_unfiltered})",
                f"Filtered\n({n_filtered})",
                f"Drought\nSignature\n({n_signature})"]
            plt.xticks(ticks=[0, 1, 2], labels=xtick_labels)

            plt.ylim(0, 1)
            plt.title(f"#{num_samples} samples")
            plt.ylabel("Accuracy")
            plt.legend(loc="lower right")
            
            # Save plot
            out_file = os.path.join(config["plotting_dir"], f"rfs_{study_name}_l{level}.svg")
            plt.savefig(out_file, dpi=600, transparent=True, bbox_inches='tight')
            plt.savefig(out_file.replace("svg", "png"), dpi=600, transparent=True, bbox_inches='tight')
            plt.close()
            print(f"[INFO] Classifier plot saved to {out_file}")

            # select top/bottom 10 features by logFC
            sorted_feats = sorted(
                diff_abundance_results.items(),
                key=lambda x: x[1]['logFC']
            )
            import re
            sorted_feats = [x for x in sorted_feats if  bool(re.fullmatch(r'[A-Za-z_;]+', x[0]))]
            keys = [k for k,_ in sorted_feats[:20]] + [k for k,_ in sorted_feats[-20:]]

            df = unfiltered_df[unfiltered_df["Height"].notnull()].copy()
            feature_intersection = list(set(df.columns) & set(keys))
            X = df[feature_intersection]

            # binary covariates
            X["Treatment"] = df["Treatment"].eq("Drought").astype(int)
            X["Inoculum"]  = df["Inoculum"].eq("Drought legacy").astype(int)
            y = df["Height"].astype(float)
            print(y)

            scalar = StandardScaler()
            enet = ElasticNetCV(cv=10, random_state=42)
            pipeline = Pipeline([('transformer', scalar),
                                 ('estimator', enet)])

            print(X)
            print(y)
            # out‑of‑fold ElasticNet predictions
            y_pred = cross_val_predict(pipeline, X, y, cv=10, n_jobs=-1)

            # spearman
            rho, pval = pearsonr(y, y_pred)

            # plot
            plt.figure(figsize=(5,5))
            plt.scatter(y_pred, y, alpha=0.7, edgecolor='k', color='black')
            mn, mx = y.min(), y.max()
            plt.plot([mn, mx], [mn, mx], 'r--')
            plt.xlabel("Predicted Height")
            plt.ylabel("Observed Height")
            plt.title(f"Correlation ρ = {rho:.2f} ({len(X)} samples)")
            plt.tight_layout()

            out = os.path.join(config["plotting_dir"], f"enet_cv_pred_{study_name}_l{level}.png")
            plt.savefig(out, dpi=600, transparent=True)
            plt.savefig(out.replace("png", "svg"), transparent=True)
            plt.close()
            print(f"[INFO] ElasticNet CV plot saved to {out}")

            colors = df["Treatment"].map({"Drought": "black"}).fillna("gray")
            plt.figure(figsize=(5,5))
            plt.scatter(y_pred, y, c=colors, alpha=0.7)
            plt.plot([mn, mx], [mn, mx], 'r--')
            plt.xlabel("Predicted Height")
            plt.ylabel("Observed Height")
            plt.title(f"Correlation ρ = {rho:.2f} ({len(X)} samples)")
            plt.tight_layout()
            plt.savefig(out.replace(".png", "_treatment.png"), dpi=600, transparent=True)
            plt.close()

            # Differential abundance analysis for inocula
            from rpy2 import robjects
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.packages import importr
            from rpy2.robjects import Formula

            # limma input requires separate metadata and transposed abundance matrix
            f = os.path.join(comp_dir, "feature-table-unfiltered.tsv")
            count_df = pd.read_csv(f, sep='\t', index_col=0)
            count_df = count_df.join(metadata_df)
            count_df = count_df[count_df['Sample type'] == 'Plant-associated']
            taxonomic_features = [col for col in count_df.columns if "p__" in col]
            metadata = count_df[["Inoculum", "Treatment"]]
            count_df = count_df[taxonomic_features]
            metadata["Inoculum"] = metadata["Inoculum"].apply(
                lambda x: "DroughtLegacy" if x == "Drought legacy" else "NotDroughtLegacy").astype("category")
            metadata["Treatment"] = metadata["Treatment"].astype("category")
            abundance_r = pandas2ri.py2rpy(count_df.T)
            metadata_r = pandas2ri.py2rpy(metadata)
            robjects.globalenv['v'] = abundance_r
            robjects.globalenv['metadata'] = metadata_r

            print(metadata)

            # Fit limma
            robjects.r('''
            library(edgeR)
            library(limma)

            dge <- DGEList(counts = v)
            dge <- calcNormFactors(dge)

            # Control is the baseline
            metadata$Inoculum <- relevel(metadata$Inoculum, ref="NotDroughtLegacy")
            design <- model.matrix(~ Treatment + Inoculum, data=metadata)
            
            # voom transformation
            v_voom <- voom(dge, design, plot=FALSE)

            # fit model
            fit <- lmFit(v_voom, design)
            fit <- eBayes(fit)

            tt <- topTable(fit, coef="InoculumDroughtLegacy", p=0.05, number=Inf, adjust.method="BH")
            ''')

            tt_df = pandas2ri.rpy2py(robjects.globalenv['tt'])

            results = []
            for feature, row in tt_df.iterrows():
                if feature in keys:
                    results.append([feature, row['logFC']])

            inoculum_diff_abundance_per_study[study_name] = results

        print(inoculum_diff_abundance_per_study)

        # 1) build DataFrame
        # --------------------------------
        # dict: study → list of [feature, logFC, og_logFC]
        # e.g. {'study1': [['d__…;p__…;…;g__X;…', 1.2], ...], …}
        data = {}
        for study, lst in inoculum_diff_abundance_per_study.items():
            data[study] = {feat: logfc for feat, logfc, og_logFC in lst}

        df_fc = pd.DataFrame.from_dict(data, orient='index').fillna(0)

        # 2) sort columns by phylum tag "p__…"
        # --------------------------------------
        def logfc_key(feat):
            return -diff_abundance_results.get(feat, {'logFC': 0})['logFC']  # negative for descending order

        sorted_cols = sorted(df_fc.columns, key=logfc_key)

        print(df_fc)

        # 3) genus labels for x‐ticks
        # --------------------------------
        def genus_label(feat):
            for part in feat.split(';'):
                if part.startswith('g__'):
                    return part.replace('g__', '')
            return feat  # fallback

        genus_labels = [genus_label(f) for f in df_fc.columns]

        # 4) compute max_abs and build cmap
        # -----------------------------------
        max_abs = np.abs(df_fc.values).max()

        cmap = LinearSegmentedColormap.from_list(
            'black_white_orange',
            ['black', 'white', '#fc8135']
        )

        # 5) plot
        # -----------------------------------
        plt.figure(figsize=(len(df_fc.columns)*0.3 + 2, len(df_fc)*0.4 + 2))
        sns.heatmap(
            df_fc,
            cmap=cmap,
            center=0,
            vmin=-max_abs,
            vmax=+max_abs,
            cbar_kws={'label': 'logFC'},
            xticklabels=genus_labels,
            yticklabels=df_fc.index
        )
        plt.xticks(rotation=90, fontsize='small')
        plt.yticks(rotation=0)
        plt.xlabel('Genus')
        plt.ylabel('Study')
        plt.tight_layout()

        out_hm = os.path.join(config["plotting_dir"], f"inoculum_logFC_heatmap_l{level}.svg")
        plt.savefig(out_hm, dpi=600, transparent=True)
        plt.close()
        print(f"[INFO] Heatmap saved to {out_hm}")
            
        
        # if level == 6:
        #     plot_bubbles(merged_file,
        #                  config["plotting_dir"])

    #     plot_num_samples(merged_file,
    #                      config["plotting_dir"])

    #     plot_taxonomy(level,
    #                   merged_file,
    #                   config["plotting_dir"])
        
    #     sparsity_heatmap(level,
    #                      level_name,
    #                      merged_file,
    #                      config["plotting_dir"])
    #     sparsity_heatmap(level,
    #                      level_name,
    #                      merged_file,
    #                      config["plotting_dir"],
    #                      log_scale=True)
    #     sparsity_heatmap(level,
    #                      level_name,
    #                      merged_file,
    #                      config["plotting_dir"],
    #                      log_scale=True,
    #                      shuffle=True)
        # plot_pcoa(level,
        #           level_name,
        #           merged_file,
        #           config)

if __name__ == "__main__":
    main()
