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

        # Differential abundance analysis
        diff_abundance_results = run_limma_diff_abundance_treatment(level,
                                                                    merged_file_counts)
        # Filter out "uncultured" bacteria
        # (this should have been done at the beginning so I need to fix this)
        filtered_results = {key: diff_abundance_results[key]
                            for key in diff_abundance_results
                            if "uncultured" not in key}
        diff_abundance_results = filtered_results

        processed_inoculum_files_counts = dict()
        metadata_inoculum_files = dict()
        diff_abundance_results_inoculum = dict()
        for study in config["inoculum_studies"]:
            study_path = os.path.join(config["data_path"],
                                      study)
            filtered_tsv, filtered_tsv_counts = preprocess_and_filter(study_path,
                                                                      study,
                                                                      level,
                                                                      config["levels"][level][0] + "__")
            processed_inoculum_files_counts[study] = filtered_tsv_counts

            # process metadata
            metadata_file = process_metadata(study_path,
                                             study,
                                             config[study])
            metadata_inoculum_files[study] = metadata_file

            merged = pd.read_csv(filtered_tsv_counts, sep='\t', index_col=0)
            merged = merged.join(pd.read_csv(metadata_file, sep='\t', index_col=0))
            merged = merged[merged['Sample type']   == 'Plant-associated']
            save_as = os.path.join(study_path,
                                   "merged.tsv")
            merged.to_csv(save_as, sep='\t')

            model = "Treatment"
            diff_abundance_results_inoculum[study] = run_limma_diff_abundance_treatment(level,
                                                                                        save_as,
                                                                                        p_val=0.05,
                                                                                        formula=model)

        plotted_labels = plot_lfc_diff_abundance(diff_abundance_results,
                                                 level,
                                                 config)

        # This method re-does the processing and should be changed
        run_external_signature_validation(diff_abundance_results,
                                          config,
                                          level)

        plot_diff_abundance_comparison(list(plotted_labels),
                                       diff_abundance_results,
                                       diff_abundance_results_inoculum)
        
        # plot_pcoa(level,
        #           level_name,
        #           merged_file,
        #           config)

if __name__ == "__main__":
    main()
