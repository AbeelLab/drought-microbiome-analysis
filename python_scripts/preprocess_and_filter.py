import pandas as pd
import subprocess
import os
import utils
import numpy as np

# Total sum scaling of features per row (sample)
def normalize_row(row,
                  taxonomic_features,
                  total_sum=100):
    total = row[taxonomic_features].sum()
    if total != 0:
        row[taxonomic_features] = row[taxonomic_features] / total
        row[taxonomic_features] *= total_sum
    return row

def filter_features(study_name,
                    study_path,
                    df,
                    level,
                    min_abundance=0.005,
                    min_prevalence=0.05):
    initial_count=len(df.columns)
    
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

    initial_count_no_unassigned = len(taxonomic_features)

    # Re-normalize
    df = df[["Study"] + to_keep]
    df = df.apply(lambda row: normalize_row(row, to_keep), axis=1)

    # save feature counts as .tsv
    save_stats_as = os.path.join(study_path,
                                 f"feature_filtering_stats_l{level}.tsv")
    f = open(save_stats_as,'w')
    f.write('study\tinitial\tafter_unassigned\tafter_prevalence_abundance\n')
    f.write(f'{study_name}\t{initial_count}\t{initial_count_no_unassigned}\t{len(to_keep)}')
    f.close()
    
    return df

# max_unassigned is from 0 to 100
# we pass a dataframe with samples normalized between 0 and 100
def filter_samples(study_name,
                   study_path,
                   df,
                   taxonomic_features,
                   level,
                   max_unassigned = 50,
                   max_eukaryotes=10):
    # Number of samples at the beginning
    initial_count = len(df)
    
    # Remove samples with too many Unassigned taxa
    unassigned_col = [col for col in taxonomic_features
                      if "Unassigned" in col]
    for col in unassigned_col:
        df = df[~(df[col] > max_unassigned)]
    after_unassigned_count = len(df)

    # Remove samples with high percentages of Eukaryotes
    # These should be removed since before downloading,
    # But for some studies such samples are not labelled
    # So it is not possible
    eukaryotic_col = [col for col in taxonomic_features
                      if "Eukaryota" in col]
    for col in eukaryotic_col:
        df = df[~(df[col] > max_eukaryotes)]
    after_eukaryotic_count = len(df)    
    
    # Remove samples with only zeros
    df = df[(df[taxonomic_features] != 0).any(axis=1)]
    final_count = len(df)

    # save counts as .csv
    save_stats_as = os.path.join(study_path,
                                 f"sample_filtering_stats_l{level}.tsv")
    f = open(save_stats_as,'w')
    f.write('study\tinitial\tafter_unassigned\tafter_eukaryotes\tafter_zeros\n')
    f.write(f'{study_name}\t{initial_count}\t{after_unassigned_count}\t{after_eukaryotic_count}\t{final_count}')
    f.close()
    
    return df

# comp_dir = directory where the extracted taxonomic composition table is stored
def biom_to_tsv(study_name,
                comp_dir):
    data_dir = utils.get_qiime_extract_dir(comp_dir)
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
    comp_dir = os.path.join(study_path, f"composition_table_l{level}")

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
    df = filter_samples(study_name,
                        study_path,
                        df,
                        taxonomic_features,
                        level)

    # Save df before feature filtering
    tsv_file_unfiltered = os.path.join(comp_dir,
                                       "feature-table-before-feature-filtering.tsv")
    df.to_csv(tsv_file_unfiltered, sep='\t')
    print(f"[INFO] Saved unfiltered feature table to: {tsv_file_unfiltered}")

    # Filter features
    df = filter_features(study_name,
                         study_path,
                         df,
                         level)

    # Save filtered df
    tsv_file_filtered = os.path.join(comp_dir, "feature-table-filtered.tsv")
    df.to_csv(tsv_file_filtered, sep='\t')
    print(f"[INFO] Saved filtered feature table to: {tsv_file_filtered}")

    return tsv_file_filtered

def process_metadata(study_path,
                     study_name,
                     study_info):
    # This already exists and should have been merged with the
    # supplemental data (if it exists) during the initial filtering
    metadata_file = os.path.join(study_path,
                                 "metadata.tsv")

    # Load metadata
    # Index column (ID) represents the sample IDs
    df = pd.read_csv(metadata_file,
                     sep='\t',
                     index_col=0)

    # rename treatment columns consistently across studies
    if "treatment_col" in study_info:
        # TO DO: this can probably be done better
        # but for some studies (simmons2020drought), we already have a treatment column
        # which gets duplicated
        if "Treatment" in df.columns:
            df.drop(columns=["Treatment"], inplace=True)
            
        df.rename(columns={study_info["treatment_col"]: "Treatment"}, inplace=True)
        df["Treatment"] = df["Treatment"].replace(study_info["treatments"],
                                                  regex=True)
    else:
        df["Treatment"] = np.nan

    # rename host columns consistently across studies
    if "host_col" in study_info:
        df["Host (original NCBI)"] = df[study_info["host_col"]]
        df["Host (specific)"] = df[study_info["host_col"]].replace(study_info["hosts_ungrouped"],
                                                                  regex=True)
        df.rename(columns={study_info["host_col"]: "Host"}, inplace=True)
        df["Host"] = df["Host"].replace(study_info["hosts"],
                                        regex=True)
    else:
        df["Host"] = np.nan
        df["Host (specific)"] = np.nan

    # rename inoculum columns consistenly across studies
    if "inoculum_col_name" in study_info:
        df.rename(columns={study_info["inoculum_col_name"]: "Inoculum"}, inplace=True)
        df["Inoculum"] = df["Inoculum"].replace({study_info["normal_inoculum"]: "Normal",
                                                 study_info["dry_inoculum"]: "Drought-legacy"},
                                                regex=True)
    else:
        df["Inoculum"] = np.nan

    df = df[["Treatment", "Host", "Host (specific)", "Host (original NCBI)", "Inoculum"]]

    metadata_file_processed = os.path.join(study_path, "processed_metadata.tsv")
    df.to_csv(metadata_file_processed, sep='\t')
    print(f"[INFO] Saved processed metadata to: {metadata_file_processed}")
    
    return metadata_file_processed
