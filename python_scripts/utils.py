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
                   NCBI_metadata,
                   supplemental_metadata,
                   link_NCBI_metadata,
                   link_supplemental_metadata,
                   save_as):
    NCBI_df = pd.read_csv(NCBI_metadata,
                          sep='\t',
                          index_col="ID")
    supplemental_df = pd.read_csv(supplemental_metadata,
                                  sep='\t')

    merged_df = NCBI_df.reset_index().merge(supplemental_df,
                                            left_on=link_NCBI_metadata,
                                            right_on=link_supplemental_metadata).set_index('ID')

    print(merged_df)

    print(merged_df.columns)
    merged_df.to_csv(save_as,
                     sep="\t",
                     index_label="ID")

    return
