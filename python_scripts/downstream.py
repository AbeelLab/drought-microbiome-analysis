import yaml
import os
import argparse

from preprocess_and_filter import preprocess_and_filter, process_metadata
from merge_datasets import merge_datasets
from plotting import plot_pca, sparsity_heatmap, plot_num_samples, plot_taxonomy

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
            feature_file = os.path.join(study_path, f"feature_filtering_stats_l{level}.tsv")
            sample_file = os.path.join(study_path, f"sample_filtering_stats_l{level}.tsv")
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
    
    for level in config["levels"]:
        print(f"[INFO] Level: {level}")
        processed_files[level] = dict()
        
        for study in config["drought_studies"]:
            print(f"[INFO] Study: {study}")
            study_path = os.path.join(config["data_path"], study)
            processed_tsv = preprocess_and_filter(study_path, study, level)
            processed_files[level][study] = processed_tsv
        
        print()
    
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
        metadata_file = process_metadata(study_path, study, config[study])
        metadata_files[study] = metadata_file

    for level, level_name in config["levels"].items():
        merge_stats()

    
    for level, level_name in config["levels"].items():
        merged_file = merge_datasets(config["data_path"], level, processed_files[level], metadata_files)

        plot_num_samples(merged_file,
                         config["plotting_dir"])

        plot_taxonomy(level,
                      merged_file,
                      config["plotting_dir"])
        
        sparsity_heatmap(level,
                         level_name,
                         merged_file,
                         config["plotting_dir"])
        sparsity_heatmap(level,
                         level_name,
                         merged_file,
                         config["plotting_dir"],
                         log_scale=True)
        sparsity_heatmap(level,
                         level_name,
                         merged_file,
                         config["plotting_dir"],
                         log_scale=True,
                         shuffle=True)
        plot_pca(level,
                 level_name,
                 merged_file,
                 config["plotting_dir"])

if __name__ == "__main__":
    main()
