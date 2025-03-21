import argparse
import os
import re
import pandas as pd
import subprocess

# processed_files is a dict
# study_id: path_to_tsv_file
def merge_datasets(data_path,
                   level,
                   processed_files):
    # load dataframes
    # indices are sample IDs
    dfs = {study: pd.read_csv(processed_files[study],
                              sep='\t',
                              index_col=0)
           for study in processed_files}
    
    # resulting columns are the union of columns from individual datasets
    merged_df = pd.concat(dfs.values(),
                          axis=0).fillna(0)

    # sort columns (for easier viewing)
    taxonomic_features = merged_df.columns.drop("Study")
    ordered_columns = ["Study"] + taxonomic_features.tolist()
    merged_df = merged_df[ordered_columns]

    merged_out_file = os.path.join(data_path, f"merged_taxonomy_normalized_l{level}.tsv")
    merged_df.to_csv(merged_out_file, sep='\t')
    print(f"[INFO] Merged normalized table saved to {merged_out_file}")

    return merged_out_file
