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
                    level_shortcut,
                    min_abundance=0.005,
                    min_prevalence=0.05):
    initial_count=len(df.columns)
    
    # Remove taxa not assigned at this level
    to_keep = [col for col in df.columns
               if (col == "Study" or level_shortcut in col)]
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
                          level,
                          level_shortcut,
                          samples_to_keep=None):
    comp_dir = os.path.join(study_path, f"composition_table_l{level}")

    df = biom_to_tsv(study_name, comp_dir)

    # Normalize taxonomic features (total sum scaling)
    # We filter out ASVs without phylum-level annotations
    # So taxonomic features should contain substring p__
    taxonomic_features = [col for col in df.columns
                          if "p__" in col]
    df = df.apply(lambda row: normalize_row(row, taxonomic_features),
                  axis=1)

    # Discard samples based on annotation quality
    if samples_to_keep is not None:
        pass
    
    # Save df before feature filtering
    tsv_file_unfiltered = os.path.join(comp_dir,
                                       "feature-table-before-feature-filtering.tsv")
    df.to_csv(tsv_file_unfiltered, sep='\t')
    print(f"[INFO] Saved unfiltered feature table to: {tsv_file_unfiltered}")

    # Filter features
    df = filter_features(study_name,
                         study_path,
                         df,
                         level,
                         level_shortcut)

    # Save filtered df
    tsv_file_filtered = os.path.join(comp_dir, "feature-table-filtered.tsv")
    df.to_csv(tsv_file_filtered, sep='\t')
    print(f"[INFO] Saved filtered feature table to: {tsv_file_filtered}")

    return tsv_file_filtered

def process_metadata(study_path,
                     study_name,
                     study_info):
    metadata_file = os.path.join(study_path, "metadata.tsv")
    df = pd.read_csv(metadata_file, sep='\t', index_col=0)

    # Map treatment columns
    if "treatment_col" in study_info:
        if "Treatment" in df.columns:
            df.drop(columns=["Treatment"], inplace=True)
        df.rename(columns={study_info["treatment_col"]: "Treatment"}, inplace=True)
        df["Treatment"] = df["Treatment"].replace(
            study_info.get("treatments", {}), regex=True
        )
    else:
        df["Treatment"] = np.nan

    # Handle drought timepoints: extract TP only for relevant stages
    if "drought_timepoints" in study_info \
       and "timepoint_col" in study_info \
       and "tp_regex" in study_info:
        df["TP"] = np.nan
        for stage, bounds in study_info["drought_timepoints"].items():
            mask_stage = df["Treatment"] == stage
            if not mask_stage.any():
                continue
            col = study_info["timepoint_col"]
            series = df.loc[mask_stage, col]
            
            # extract timepoints with regex
            tp_extracted = series.str.extract(study_info["tp_regex"])
            tp_values = tp_extracted["tp"]
            df.loc[mask_stage, "TP"] = pd.to_numeric(tp_values, errors="coerce")

            # Determine drought vs control based on TP bounds
            if "max_tp" in bounds:
                drought_mask = mask_stage & (df["TP"] <= bounds["max_tp"])
            elif "min_tp" in bounds:
                drought_mask = mask_stage & (df["TP"] >= bounds["min_tp"])
            else:
                continue
            control_mask = mask_stage & ~drought_mask
            df.loc[drought_mask, "Treatment"] = "Drought"
            df.loc[control_mask, "Treatment"] = "Control"
        # Clean up TP column
        df.drop(columns=["TP"], inplace=True)

    # Map host columns
    if "host_col" in study_info:
        df["Host (original NCBI)"] = df[study_info["host_col"]]
        df["Host (specific)"] = (
            df[study_info["host_col"]]
              .replace(study_info.get("hosts_ungrouped", {}), regex=True)
        )
        df.rename(columns={study_info["host_col"]: "Host"}, inplace=True)
        df["Host"] = df["Host"].replace(
            study_info.get("hosts", {}), regex=True
        )
    else:
        df["Host"] = np.nan
        df["Host (specific)"] = np.nan

    # Map inoculum columns
    if "inoculum_col_name" in study_info:
        df.rename(columns={study_info["inoculum_col_name"]: "Inoculum"}, inplace=True)
        df["Inoculum"] = df["Inoculum"].replace({
            study_info.get("normal_inoculum"): "Normal",
            study_info.get("dry_inoculum"): "Drought-legacy"
        }, regex=True)
    else:
        df["Inoculum"] = np.nan

    # Add study-specific metadata columns
    df["Study (full name)"] = study_info.get("full_name", np.nan)
    df["Location"] = study_info.get("location", np.nan)
    df["Soil type"] = study_info.get("soil_type", np.nan)
    df["Primers"] = study_info.get("primers", np.nan)
    df["Regions"] = study_info.get("regions", np.nan)

    # Select and order final columns
    final_cols = [
        "Study (full name)", "Treatment", "Host", "Host (specific)", "Host (original NCBI)",
        "Inoculum", "Location", "Soil type", "Primers", "Regions"
    ]
    df = df[final_cols]

    # Save processed metadata
    processed_file = os.path.join(study_path, "processed_metadata.tsv")
    df.to_csv(processed_file, sep='\t')
    print(f"[INFO] Saved processed metadata to: {processed_file}")
    return processed_file
