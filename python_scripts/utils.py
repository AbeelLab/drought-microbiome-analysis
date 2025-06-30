import numpy as np
import os
import pandas as pd
import subprocess

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

# comp_dir = directory where the extracted taxonomic composition table is stored
def biom_to_tsv(comp_dir):
    data_dir = get_qiime_extract_dir(comp_dir)
    biom_file = os.path.join(data_dir, "feature-table.biom")

    tsv_file_unfiltered = os.path.join(comp_dir,
                                       "feature-table-unfiltered.tsv")

    # Need to activate biom environment...
    subprocess.run(["conda", "run", "-n", "biom", "biom", "convert",
                    "-i", biom_file, "-o", tsv_file_unfiltered, "--to-tsv"], check=True)

    # Remove the first line
    # biom convert adds a useless line at the beginning of the file
    with open(tsv_file_unfiltered, 'r') as f:
        lines = f.readlines()[1:]
    with open(tsv_file_unfiltered, 'w') as f:
        f.writelines(lines)

    # Transpose so rows are samples and columns are features (taxa)
    df = pd.read_csv(tsv_file_unfiltered,
                     sep='\t',
                     index_col=0).T

    df.to_csv(tsv_file_unfiltered,
              sep='\t')
    print(f"[INFO] Saved unfiltered feature table to: {tsv_file_unfiltered}")

    return df

def merge_metadata(study, config, NCBI_metadata, supplemental_metadata, save_as):
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

def process_metadata(study_path,
                     study_name,
                     study_info):
    metadata_file = os.path.join(study_path, "metadata.tsv")
    df = pd.read_csv(metadata_file, sep='\t', index_col=0)

    # Map treatment columns
    if "treatment_col" in study_info:
        if "Treatment" in df.columns:
            df.drop(columns=["Treatment"], inplace=True)
        # preserve original column
        df['Treatment'] = df[study_info['treatment_col']]
        df['Treatment'] = df['Treatment'].replace(
            study_info.get('treatments', {}), regex=True)
    else:
        df['Treatment'] = np.nan

    # Handle drought timepoints: extract TP only for relevant stages
    if "drought_timepoints" in study_info \
       and "timepoint_col" in study_info \
       and "tp_regex" in study_info:
        df['TP'] = np.nan
        for stage, bounds in study_info['drought_timepoints'].items():
            mask_stage = df['Treatment'] == stage
            if not mask_stage.any():
                continue
            col = study_info['timepoint_col']
            series = df.loc[mask_stage, col]
            
            # extract timepoints with regex
            tp_extracted = series.str.extract(study_info['tp_regex'])
            tp_values = tp_extracted['tp']
            df.loc[mask_stage, 'TP'] = pd.to_numeric(tp_values, errors='coerce')

            # Determine drought vs control based on TP bounds
            if 'max_tp' in bounds:
                drought_mask = mask_stage & (df['TP'] <= bounds['max_tp'])
            elif 'min_tp' in bounds:
                drought_mask = mask_stage & (df['TP'] >= bounds['min_tp'])
            else:
                continue
            control_mask = mask_stage & ~drought_mask
            df.loc[drought_mask, 'Treatment'] = 'Drought'
            df.loc[control_mask, 'Treatment'] = 'Control'
        # Clean up TP column
        df.drop(columns=['TP'], inplace=True)

    # Map host columns
    if "host_col" in study_info:
        df['HostSpecific'] = df[study_info['host_col']]
        df['HostSpecific'] = df['HostSpecific'].replace(
            study_info.get('hosts_ungrouped', {}), regex=True)
        df['Host'] = df[study_info['host_col']]
        df['Host'] = df['Host'].replace(
            study_info.get('hosts', {}), regex=True
        )
    else:
        df['Host'] = np.nan
        df['HostSpecific'] = np.nan

    # Map root compartment columns
    if "root_compartment_col" in study_info:
        df['RootCompartment'] = df[study_info['root_compartment_col']]
        df['RootCompartment'] = df['RootCompartment'].replace(
            study_info.get('root_compartments', {}), regex=True)
    else:
        df['RootCompartment'] = np.nan    

    # Remove NaN hosts (those not mapped)
    df = df[df['HostSpecific'] != 'None']

    # Map inoculum columns
    if "inoculum_col" in study_info:
        df['Inoculum'] = df[study_info['inoculum_col']]
        df['Inoculum'] = df['Inoculum'].replace(
            study_info.get('inocula', {}), regex=True)

        df['InoculumSubtypes'] = df[study_info['inoculum_col']]
        df['InoculumSubtypes'] = df['InoculumSubtypes'].replace(
            study_info.get('inocula_subtypes', {}), regex=True)
    else:
        df['Inoculum'] = np.nan
        df['InoculumSubtypes'] = np.nan

    # Sample types
    if "sample_type_col" in study_info:
        df['Sample type'] = df[study_info['sample_type_col']]
        df['Sample type'] = df['Sample type'].replace(
            study_info.get('sample_types', {}), regex=True
        )
    else:
        df['Sample type'] = np.nan

    # Add study-specific metadata columns
    df['Study'] = study_name
    df['Study (full name)'] = study_info.get('full_name', np.nan)
    df['Location'] = study_info.get('location', np.nan)
    df['SoilType'] = study_info.get('soil_type', np.nan)
    df['Primers'] = study_info.get('primers', np.nan)
    df['Regions'] = study_info.get('regions', np.nan)

    # Select and order final columns
    final_cols = [
        'Study',
        'Study (full name)',
        'Treatment',
        'Host',
        'HostSpecific',
        'Inoculum',
        'InoculumSubtypes',
        'Location',
        'SoilType',
        'Primers',
        'Regions',
        'Sample type'
    ]
    df = df[final_cols]

    return df
