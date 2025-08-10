import os
import argparse
import utils
import yaml
import utils
import pandas as pd
import shutil

def apply_filters(df, filters):
    for filter_expr in filters:
        print(filter_expr)
        filter_func = eval(filter_expr)
        df = df[df.apply(filter_func, axis=1)]
    return df

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study")
    args = parser.parse_args()

    study = args.study

    with open('config.yml') as f:
        config = yaml.safe_load(f)

    study_dir = os.path.join(config["data_path"],
                             study)
    metadata_dir = os.path.join(study_dir,
                                "metadata")
    inner_dir = utils.get_qiime_extract_dir(metadata_dir)
    metadata_file = os.path.join(inner_dir,
                                 "sra-metadata.tsv")

    # Count number of lines, minus heading
    num_samples = sum(1 for _ in open(metadata_file)) - 1
    # Save #samples to .pkl file
    utils.log_statistics(study, "Total NCBI runs", num_samples)

    # Copy file to main directory as well
    copy_to = os.path.join(study_dir,
                           "metadata.tsv")
    shutil.copyfile(metadata_file, copy_to)

    # If needed, merge with supplemental metadata
    # Overwrites metadata.tsv
    supplemental_metadata = os.path.join(study_dir,
                                         "supplemental_metadata.tsv")
    if os.path.exists(supplemental_metadata):
        utils.merge_metadata(study,
                             config,
                             metadata_file,
                             supplemental_metadata,
                             copy_to)

    df = pd.read_csv(copy_to,
                     sep="\t",
                     index_col="ID")

    utils.log_statistics(study,
                         "NCBI runs after merging with study supplemental data",
                         len(df))

    print("Available samples: ", len(df))

    # Apply individual study filters
    if study in config and "filters" in config[study]:
        filtered_df = apply_filters(df,
                                    config[study].get("filters",
                                                      []))
    else:
        print(f"No filters for study {study}")
        filtered_df = df

    # Log #samples after initial filtering
    utils.log_statistics(study,
                         "NCBI runs after initial filtering",
                         len(filtered_df))
    utils.log_statistics(study,
                         "NCBI BioSamples after initial filtering",
                         filtered_df["Biosample ID"].nunique())    
        
    filtered_df.to_csv(os.path.join(study_dir, "filtered_df.tsv"),
                       sep='\t')

    print("Available runs after filtering: ",
          len(filtered_df))
    print("Available samples after filtering: ",
          filtered_df["Biosample ID"].nunique())

    # Merge runs belonging to one sample by taking one (the one with max number of bases)
    filtered_df.sort_values('Bases', ascending=False, inplace=True)
    filtered_df.drop_duplicates(subset='Biosample ID', keep="first", inplace=True)

    assert len(filtered_df) == filtered_df["Biosample ID"].nunique()
    
    # Save only IDs (one column) to filtered_sras.tsv
    # Rename ID to id
    ids_df = filtered_df.reset_index()[["ID"]].rename(columns={"ID": "id"})
    ids_df.to_csv(os.path.join(study_dir, "filtered_sras.tsv"),
                  sep='\t',
                  index=False)
    

if __name__ == "__main__":
    main()
