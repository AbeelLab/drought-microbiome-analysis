import yaml
import os
import argparse

from preprocess_and_filter import preprocess_and_filter, process_metadata
from merge_datasets import merge_datasets
from plotting import plot_pca, sparsity_heatmap

with open('config.yml') as f:
        config = yaml.safe_load(f)

def run_preprocessing():
    processed_files = dict()
    
    for level in config["levels"]:
        print(f"[INFO] Level: {level}")

        processed_files[level] = dict()
        
        for study in config["studies"]:
            print(f"[INFO] Study: {study}")
            study_path = os.path.join(config["data_path"],
                                      study)

            processed_tsv = preprocess_and_filter(study_path,
                                                  study,
                                                  level)
            processed_files[level][study] = processed_tsv
        
        print()

    with open(config["processed_files"], 'w') as f:
        yaml.dump(processed_files,
                  f,
                  default_flow_style=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--preprocess", action=argparse.BooleanOptionalAction)
    args = parser.parse_args()

    preprocess = args.preprocess

    # To avoid waiting for preprocessing each time
    if (not os.path.exists(config["processed_files"])) or preprocess:
        run_preprocessing()

    with open(config["processed_files"]) as f:
        processed_files = yaml.safe_load(f)

    metadata_dfs = dict()
    # Process metadata per study
    for study in config["studies"]:
        if config["studies"][study] == "No data":
                continue
        
        study_path = os.path.join(config["data_path"],
                                  study)
        metadata_df = process_metadata(study_path,
                                       study,
                                       config["studies"][study])
        metadata_dfs[study]=metadata_df
        

    for level, level_name in config["levels"].items():
        # Rarefaction analysis
        # plot_rarefaction(processed_files[level])
        # what is the cap?
        # (i.e. how many taxa, per level, are documented in SILVA?)
        
        # Merge per level
        merged_file = merge_datasets(config["data_path"],
                                     level,
                                     processed_files[level])

        # Plot data heatmap
        sparsity_heatmap(level,
                         level_name,
                         merged_file,
                         config["plotting_dir"])
        # Log-scale colors
        sparsity_heatmap(level,
                         level_name,
                         merged_file,
                         config["plotting_dir"],
                         log_scale=True)
        # Log-scale colors and column shuffling
        sparsity_heatmap(level,
                         level_name,
                         merged_file,
                         config["plotting_dir"],
                         log_scale=True,
                         shuffle=True)
        
        # Plot PCA
        plot_pca(level,
                 level_name,
                 merged_file,
                 config["plotting_dir"])
    

if __name__ == "__main__":
    main()
