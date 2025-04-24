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
    # Remove incomplete records, or otherwise records not mapped
    # to a group of Hosts
    merged_df = merged_df.dropna(subset=['Host'])

    # sort columns (for easier viewing)
    metadata = [col for col in merged_df.columns
                if "p__" not in col]
    taxonomic_features = [col for col in merged_df.columns
                          if "p__" in col]

    # Save separately as well
    # technically these are already both saved separately, but
    # as a sanity check should be saved after the join (in case some
    # samples without metadata should get dropped)
    merged_metadata_file = os.path.join(data_path, f"merged_metadata_l{level}.tsv")
    merged_df[metadata].to_csv(merged_metadata_file, sep='\t', index_label="ID")
    
    merged_taxonomy_file = os.path.join(data_path, f"merged_only_taxa_l{level}.tsv")
    # Divide by 100 because QIIME2 RelativeFrequency artifact requires values to sum up to 1
    # Then we should transpose and rename the index #OTU ID for biom convert to work?
    taxa = merged_df[taxonomic_features].div(100).T
    taxa.index.name = '#OTU ID'
    with open(merged_taxonomy_file, 'w') as f:
        f.write('# Constructed from biom file\n')
        taxa.to_csv(f, sep='\t')
    
    ordered_columns = metadata + taxonomic_features
    merged_df = merged_df[ordered_columns]

    merged_out_file = os.path.join(data_path, f"merged_taxonomy_normalized_l{level}.tsv")
    merged_df.to_csv(merged_out_file, sep='\t')
    print(f"[INFO] Merged normalized table saved to {merged_out_file}")

    return merged_out_file
