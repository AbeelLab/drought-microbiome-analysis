import numpy as np
import pandas as pd

# Column names
treatment_col = "Treatment"
host_col = "HostSpecific"
grouped_host_col = "Host"
root_compartment_col = "RootCompartment"
inoculum_col = "Inoculum"
inoculum_subtypes_col = "InoculumSubtype"
is_plant_associated_col = "IsPlantAssociated"
study_col = "StudyID"
full_study_col = "StudyName"
primer_col = "Primers"
regions_col = "Regions"

def common_processing(study_path,
                      study_name,
                      study_info):
    metadata_file = os.path.join(study_path, "metadata.tsv")
    df = pd.read_csv(metadata_file, sep='\t', index_col=0)

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
        
    # Remove NaN hosts (those not mapped)
    # This is because some hosts appear in very small groups
    df = df[df[host_col] != 'None']

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
    df[region_col] = study_info.get('regions', np.nan)

    # CALL PER-DATASET PROCESSING FUNCTION HERE
    if not None:
        df = process_x(...)

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
        region_col,
        is_plant_associated_col
    ]
    df = df[final_cols]

    return df

def process_xu2018drought(df):
    # First field experiment
    def extract_tp(tp_string):
        tp_component = [x for x in tp_string.split("_")
                        if "TP" in x][0]
        return int(tp_component[2:]) * 7

    for i, row in df.iterrows():
        if row["Treatment"] == "Pre_flowering":
            # Remove timepoints before 2nd week (seedling development)
            if extract_tp(row["Title"]) <= 2:
                df.at[i, "Treatment"] = np.nan
            elif extract_tp(row["Title"]) <= 8:
                df.at[i, "Treatment"] = "Drought"
            else:
                df.at[i, "Treatment"] = "Control"
        if row["Treatment"] == "Post_flowering":
            # Remove timepoints before 2nd week (seedling development)
            if extract_tp(row["Title"]) <= 2:
                df.at[i, "Treatment"] = np.nan
            elif extract_tp(row["Title"]) >= 10:
                #  Flowering happens at week 9, and post-flowering drought starts in week 10
                df.at[i, "Treatment"] = "Drought"
            else:
                df.at[i, "Treatment"] = "Control"

    # Remove NaNs (seedling development)
    df = df[df["Treatment"].notna()]

    return df
