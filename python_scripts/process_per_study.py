import importlib
import numpy as np
import pandas as pd
import os
import math

# Column names
treatment_col = "Treatment"
host_col = "HostSpecific"
grouped_host_col = "Host"
root_compartment_col = "RootCompartment"
inoculum_col = "Inoculum"
inoculum_subtype_col = "InoculumSubtype"
is_plant_associated_col = "IsPlantAssociated"
study_col = "StudyID"
full_study_col = "StudyName"
primer_col = "Primers"

def load_function(full_name):
    module_name, func_name = full_name.rsplit('.', 1)
    module = importlib.import_module(module_name)
    return getattr(module, func_name)

def common_processing(study_path,
                      study_name,
                      study_info,
                      custom_processing=None,
                      kept_samples=None):        
    metadata_file = os.path.join(study_path, "metadata.tsv")
    df = pd.read_csv(metadata_file, sep='\t', index_col=0)

    if kept_samples is not None:
        df = df.loc[kept_samples]
        if study_name == "azarbad2022response":
            df.to_csv("godhelpme.tsv", sep="\t")

    # Map treatment columns
    if "treatment_col" in study_info:
        df[treatment_col] = df[study_info['treatment_col']]
        df[treatment_col] = df[treatment_col].replace(
            study_info.get('treatments', {}), regex=True)
    else:
        df[treatment_col] = np.nan

    # Map host columns
    if "host_col" in study_info:
        # Ungrouped hosts
        df[host_col] = df[study_info['host_col']]
        df[host_col] = df[host_col].replace(
            study_info.get('hosts_ungrouped', {}), regex=True)

        # Grouped hosts
        df[grouped_host_col] = df[study_info['host_col']]
        df[grouped_host_col] = df[grouped_host_col].replace(
            study_info.get('hosts', {}), regex=True
        )
    else:
        df[grouped_host_col] = np.nan
        df[host_col] = np.nan

    # Map root compartment columns
    if "root_compartment_col" in study_info:
        df[root_compartment_col] = df[study_info['root_compartment_col']]
        df[root_compartment_col] = df[root_compartment_col].replace(
            study_info.get('root_compartments', {}), regex=True)
    else:
        df[root_compartment_col] = np.nan    

    # Map inoculum columns
    if "inoculum_col" in study_info:
        df[inoculum_col] = df[study_info['inoculum_col']]
        df[inoculum_col] = df[inoculum_col].replace(
            study_info.get('inocula', {}), regex=True)

        df[inoculum_subtype_col] = df[study_info['inoculum_col']]
        df[inoculum_subtype_col] = df[inoculum_subtype_col].replace(
            study_info.get('inocula_subtypes', {}), regex=True)
    else:
        df[inoculum_col] = "Not available for this type of experiment"
        df[inoculum_subtype_col] = np.nan

    # Sample types
    if "sample_type_col" in study_info:
        df[is_plant_associated_col] = df[study_info['sample_type_col']]
        df[is_plant_associated_col] = df[is_plant_associated_col].replace(
            study_info.get('sample_types', {}), regex=True
        )
    else:
        df[is_plant_associated_col] = np.nan

    # Add constant study-specific metadata columns
    df[study_col] = study_name
    df[full_study_col] = study_info.get('full_name', np.nan)
    df[primer_col] = study_info.get('primers', np.nan)

    # Call custom processing pe dataset
    if custom_processing is not None:
        fn = load_function(custom_processing)
        df = fn(df, study_info)

    # Replace all hosts where root_compartment_col is "Bulk soil" with Soil
    df.loc[df[root_compartment_col] == "Bulk soil", host_col] = "Soil"
    df.loc[df[root_compartment_col] == "Bulk soil", grouped_host_col] = "Soil"

    # Select and order final columns
    final_cols = [
        study_col,
        full_study_col,
        treatment_col,
        grouped_host_col,
        host_col,
        root_compartment_col,
        inoculum_col,
        inoculum_subtype_col,
        primer_col,
        is_plant_associated_col
    ]
    df = df[final_cols]

    # Remove NaN hosts (those not mapped)
    # This is because some hosts appear in very small groups
    df = df[df[host_col] != 'None']

    return df

def process_xu2018drought(df, study_info):
    # First field experiment
    def extract_tp(tp_string):
        tp_component = [x for x in tp_string.split("_")
                        if "TP" in x][0]
        return int(tp_component[2:]) * 7

    for i, row in df.iterrows():
        if row[treatment_col] == "Pre_flowering":
            # Remove timepoints before 2nd week (seedling development)
            if extract_tp(row["Title"]) <= 2:
                df.at[i, treatment_col] = np.nan
            elif extract_tp(row["Title"]) <= 8:
                df.at[i, treatment_col] = "Drought"
            else:
                df.at[i, treatment_col] = "Control"
        if row[treatment_col] == "Post_flowering":
            # Remove timepoints before 2nd week (seedling development)
            if extract_tp(row["Title"]) <= 2:
                df.at[i, treatment_col] = np.nan
            elif extract_tp(row["Title"]) >= 10:
                #  Flowering happens at week 9,
                # and post-flowering drought starts in week 10
                df.at[i, treatment_col] = "Drought"
            else:
                df.at[i, treatment_col] = "Control"

    # Remove NaNs (seedling development)
    df = df[df[treatment_col].notna()]

    return df

def process_naylor2017drought(df, study_info):
    # Separate processing for Kearney experiments
    for i, row in df.iterrows():
        if "Kearney" in row["Geo Loc Name [sample]"]:
            if row["Name"].split("_")[2] == "D":
                df.at[i, treatment_col] = "Drought"
            elif row["Name"].split("_")[2] == "W":
                df.at[i, treatment_col] = "Control"

            if "RZ" in row["Name"]:
                df.at[i, root_compartment_col] = "Rhizosphere"
            elif "Root" in row["Name"]:
                df.at[i, root_compartment_col] = "Endosphere"
            elif "Soil" in row["Name"]:
                df.at[i, root_compartment_col] = "Bulk soil"
    for i, row in df.iterrows():
        if pd.isnull(row["Host"]):
            col = study_info["host_col_separate_experiment"]
            host_mapping = study_info["hosts_ungrouped"]
            grouped_host_mapping = study_info["hosts"]
            df.at[i, host_col] = host_mapping[df.at[i, col]]
            df.at[i, grouped_host_col] = grouped_host_mapping[df.at[i, col]]
            
    return df

def process_simmons2020drought(df, study_info):
    # Notes from supplemental metadata:
    # Sample names swapped. Fasta file reads SO-Wk9C-R3-2-3, but should be SO-Wk9D-R3-3-3.
    # Sample names swapped. Fasta file reads SO-Wk9D-R3-3-3, but should be SO-Wk9C-R3-2-3.
    
    # So we should swap the accessions for these two sample IDs
    tmp = df.loc["SRR11143478"].copy()
    df.loc["SRR11143478"] = df.loc["SRR11143463"]
    df.loc["SRR11143463"] = tmp
    return df
