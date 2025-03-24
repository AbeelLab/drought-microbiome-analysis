import os
import argparse
import utils
import yaml
import utils
import pandas as pd

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

    # Apply individual study filters
    df = pd.read_csv(metadata_file,
                     sep="\t",
                     index_col="ID")

    print("Available samples: ",
          len(df))

    filtered_df = apply_filters(df,
                                config[study].get("filters",
                                                  []))
        
    filtered_df.to_csv(os.path.join(study_dir, "filtered_df.tsv"),
                       sep='\t')

    print("Available samples after filtering: ",
          len(filtered_df))

    # Save only IDs (one column) to filtered_sras.tsv
    # Rename ID to id
    ids_df = filtered_df.reset_index()[["ID"]].rename(columns={"ID": "id"})
    ids_df.to_csv(os.path.join(study_dir, "filtered_sras.tsv"),
                  sep='\t',
                  index=False)
    

if __name__ == "__main__":
    main()
