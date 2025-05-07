import os
import pandas as pd

def get_qiime_extract_dir(parent_dir):
    # extracted file is stored in a nested directory created by qiime extract
    # (qiime export would result in a cleaner directory structure,
    # but right now it doesn't work due to permission errors on the cluster)
    subdirs = [d for d in os.listdir(parent_dir) if os.path.isdir(os.path.join(parent_dir, d))]
    # Get latest extracted directory if there are more
    nested_id_dir = max(subdirs, key=lambda d: os.path.getctime(os.path.join(parent_dir, d)))
    nested_id_dir = os.path.join(parent_dir, nested_id_dir)
    data_dir = os.path.join(nested_id_dir, "data")

    return data_dir

def merge_metadata(study,
                   config,
                   NCBI_metadata,
                   supplemental_metadata,
                   save_as):
    link_NCBI_metadata = config[study]["link_NCBI_metadata"]
    link_supplemental_metadata = config[study]["link_supplemental_metadata"]

    NCBI_df = pd.read_csv(NCBI_metadata,
                          sep='\t',
                          index_col="ID")
    supplemental_df = pd.read_csv(supplemental_metadata,
                                  sep='\t')

    # When there isn't exact matching
    if "match_regex" in config[study]:
        match_regex = config[study]["match_regex"]
        NCBI_df["match_key"] = NCBI_df[link_NCBI_metadata].str.extract(match_regex, expand=False)
        link_NCBI_metadata = "match_key"

    # Merging won't work if the columns are not the same type
    NCBI_df[link_NCBI_metadata] = NCBI_df[link_NCBI_metadata].astype(str)
    supplemental_df[link_supplemental_metadata] = supplemental_df[link_supplemental_metadata].astype(str)

    merged_df = NCBI_df.reset_index().merge(supplemental_df,
                                            how="left",
                                            left_on=link_NCBI_metadata,
                                            right_on=link_supplemental_metadata).set_index('ID')

    merged_df.to_csv(save_as,
                     sep="\t",
                     index_label="ID")

    return
