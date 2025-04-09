import argparse
import os
import re
import pandas as pd
import subprocess

# processed_files is a dict
# study_id: path_to_tsv_file
def merge_datasets(data_path,
                   level,
                   processed_files,
                   metadata_files):
    # load dataframes
    # and merge with metadata
    # indices are sample IDs
    dfs = {study: pd.read_csv(processed_files[study],
                              sep='\t',
                              index_col=0)
           for study in processed_files}

    metadata_dfs = {study: pd.read_csv(metadata_files[study],
                                       sep='\t',
                                       index_col=0)
                    for study in metadata_files}
    
    # resulting columns are the union of columns from individual datasets
    merged_df = pd.concat(dfs.values(),
                          axis=0).fillna(0)
    merged_metadata_df = pd.concat(metadata_dfs.values(),
                                   axis=0)

    merged_df = merged_df.join(merged_metadata_df)

    # sort columns (for easier viewing)
    metadata = ["Study", "Treatment", "Host", "Host (specific)", "Inoculum"]
    taxonomic_features = merged_df.columns.drop(metadata)
    ordered_columns = metadata + taxonomic_features.tolist()
    merged_df = merged_df[ordered_columns]

    merged_out_file = os.path.join(data_path, f"merged_taxonomy_normalized_l{level}.tsv")
    merged_df.to_csv(merged_out_file, sep='\t')
    print(f"[INFO] Merged normalized table saved to {merged_out_file}")

    return merged_out_file
