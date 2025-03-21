import pandas as pd
import subprocess
import os
import utils

# Total sum scaling of features per row (sample)
def normalize_row(row,
                  taxonomic_features,
                  total_sum=100):
    total = row[taxonomic_features].sum()
    if total != 0:
        row[taxonomic_features] = row[taxonomic_features] / total
        row[taxonomic_features] *= total_sum
    return row

def filter_features(df,
                    min_abundance=0.005,
                    min_prevalence=0.1):
    # Remove unassigned and uncultured taxa
    to_keep = [col for col in df.columns
               if "Unassigned" not in col and "uncultured" not in col]
    df = df[to_keep]

    # Re-normalize to enable abundance/prevalence filtering
    taxonomic_features = [col for col in df.columns if col != "Study"]
    df = df.apply(lambda row: normalize_row(row, taxonomic_features), axis=1)

    min_sample_count = int(min_prevalence * len(df))
    to_keep = [col for col in taxonomic_features
               if (df[col] >= min_abundance).sum() >= min_sample_count]

    initial_count = len(taxonomic_features)
    removed_count = initial_count - len(to_keep)
    print(f"[INFO] Removed {removed_count} features out of {initial_count}")
    print(f"[INFO] Features left: {len(to_keep)}")

    # Re-normalize
    df = df[["Study"] + to_keep]
    df = df.apply(lambda row: normalize_row(row, to_keep), axis=1)
    
    return df

def filter_samples(df,
                   taxonomic_features, 
                   max_unassigned = 50):
    # Number of samples at the beginning
    initial_count = len(df)
    
    # Remove samples with too many Unassigned taxa
    unassigned_col = [col for col in taxonomic_features
                      if "Unassigned" in col]
    if len(unassigned_col) > 0:
        unassigned_col = unassigned_col[0]
        df = df[~(df[unassigned_col] > max_unassigned)]
    after_unassigned_count = len(df)
    
    # Remove samples with only zeros
    df = df[(df[taxonomic_features] != 0).any(axis=1)]
    final_count = len(df)

    # Logging
    removed_unassigned = initial_count - after_unassigned_count
    removed_zeros = after_unassigned_count - final_count
    total_removed = removed_unassigned + removed_zeros
    
    print(f"[INFO] Initial #samples: {initial_count}")
    print(f"[INFO] Removed due to Unassigned > {max_unassigned}: {removed_unassigned}")
    print(f"[INFO] Removed empty samples: {removed_zeros}")
    print(f"[INFO] Final samples: {final_count} (Total removed: {total_removed})")
    
    return df

# comp_dir = directory where the extracted taxonomic composition table is stored
def biom_to_tsv(study_name,
                comp_dir):
    data_dir = get_qiime_extract_dir(comp_dir)
    biom_file = os.path.join(data_dir, "feature-table.biom")

    tsv_file_unfiltered = os.path.join(comp_dir, "feature-table-unfiltered.tsv")
    # Need to activate biom environment
    subprocess.run(["conda", "run", "-n", "biom", "biom", "convert",
                    "-i", biom_file, "-o", tsv_file_unfiltered, "--to-tsv"], check=True)

    # Remove the first line
    # biom convert adds a useless line at the beginning of the file
    with open(tsv_file_unfiltered, 'r') as f:
        lines = f.readlines()[1:]
    with open(tsv_file_unfiltered, 'w') as f:
        f.writelines(lines)

    # Transpose so rows are samples and columns are features (taxa)
    df = pd.read_csv(tsv_file_unfiltered, sep='\t', index_col=0).T

    # Add study name column
    df["Study"] = study_name

    df.to_csv(tsv_file_unfiltered, sep='\t')
    print(f"[INFO] Saved unfiltered feature table to: {tsv_file_unfiltered}")

    return df

# study_path = relative path to the study directory
# study_name = name of study directory
# level = taxonomic depth (from 2 (phylum) to 6 (genus))
def preprocess_and_filter(study_path,
                          study_name,
                          level):
    qiime_dir = os.path.join(study_path, "qiime-dir")
    comp_dir = os.path.join(qiime_dir, f"composition_table_l{level}")

    df = biom_to_tsv(study_name, comp_dir)

    # For now we focus on bacteria and archaea (no fungi)
    # Filter out these columns before normalizing
    # It doesn't make sense to normalize when a sample has both
    # bacterial and fungal data.
    to_remove = [col for col in df.columns
                 if ("Chloroplast" in col)
                 or ("Mitochondria" in col)
                 or ("Eukaryota" in col)]
    df.drop(columns=to_remove, inplace=True)

    # Normalize taxonomic features (total sum scaling)
    taxonomic_features = [col for col in df.columns
                          if col != "Study"]
    df = df.apply(lambda row: normalize_row(row, taxonomic_features),
                  axis=1)

    # Filter samples 
    df = filter_samples(df, taxonomic_features)

    # Filter features
    df = filter_features(df)

    # Save filtered df
    tsv_file_filtered = os.path.join(comp_dir, "feature-table-filtered.tsv")
    df.to_csv(tsv_file_filtered, sep='\t')
    print(f"[INFO] Saved filtered feature table to: {tsv_file_filtered}")

    return tsv_file_filtered

def process_metadata(study_path,
                     study_name,
                     study_info):
    metadata_dir = os.path.join(study_path, "qiime-dir", "metadata")
    metadata_dir = get_qiime_extract_dir(metadata_dir)

    metadata_file = os.path.join(metadata_dir, "sra-metadata.tsv")

    # Load metadata
    # Index column (ID) represents the sample IDs
    df = pd.read_csv(metadata_file, sep='\t', index_col=0)

    # rename column study_info["inoculum_col_name"] to "Inoculum"
    df.rename(columns={study_info["inoculum_col_name"]: "Inoculum"}, inplace=True)

    # for values in this column,
    # rename study_info["normal_inoculum"] to normal
    # rename study_info["dry_inoculum"] to drought-legacy
    df["Inoculum"] = df["Inoculum"].replace({study_info["normal_inoculum"]: "normal",
                                             study_info["dry_inoculum"]: "drought-legacy"})

    # filter such that only the inoculum column is included in the dataframe
    # and the sample IDs (the index)
    df = df[["Inoculum"]]

    metadata_file_processed = os.path.join(study_path, "processed_metadata.tsv")
    df.to_csv(metadata_file_processed, sep='\t')
    print(f"[INFO] Saved processed metadata to: {metadata_file_processed}")
    
    return df
