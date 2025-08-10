import numpy as np
import os
import pandas as pd
import pickle as pkl
import subprocess

from Bio import Phylo
from collections import defaultdict

def log_statistics(study,
                   filtering_step,
                   to_log,
                   log_file="../data/samples_log.pkl"):
    log = None
    if os.path.exists(log_file):
        log = pkl.load(open(log_file,
                            "rb"))
        if study not in log:
            log[study] = dict()
        log[study][filtering_step] = to_log
    else:
        log = dict()
        log[study] = dict()
        log[study][filtering_step] = to_log
    pkl.dump(log,
             open(log_file, "wb"))

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
                                            how="inner",
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

def log_features_per_batch(ds,
                           compartment,
                           batch_column="StudyID"):
    all_features = ds.get_counts_features_columns()
    sample_to_batch = ds.metadata_df[batch_column]
    for batch_id, sample_ids in sample_to_batch.groupby(sample_to_batch).groups.items():
        batch_df = ds.taxonomy_counts_df.loc[sample_ids, all_features]
        present_features = batch_df.columns[(batch_df > 0).any(axis=0)]
        log_statistics(f"{batch_id}",
                       f"After BE correction: {compartment}",
                       int(len(present_features)),
                       "../data/features_log.pkl")


def get_core_microbiome(ds,
                        batch_column="StudyID"):
    all_features = ds.get_counts_features_columns()
    core = set()
    sample_to_batch = ds.metadata_df[batch_column]
    for batch_id, sample_ids in sample_to_batch.groupby(sample_to_batch).groups.items():
        batch_df = ds.taxonomy_counts_df.loc[sample_ids, all_features]
        present_features = batch_df.columns[(batch_df > 0).any(axis=0)]
        core.update(set(present_features))

    return core


def compute_host_phylogenetic_distances(tree_file="../data/plant_genes/plants.nwk"):
    tree = Phylo.read(tree_file, "newick")
    terminals = tree.get_terminals()
    distance_dict = defaultdict(dict)
    for t1 in terminals:
        for t2 in terminals:
            if t1 != t2:
                dist = tree.distance(t1, t2)
                distance_dict[t1.name.replace("_", " ")][t2.name.replace("_", " ")] = dist / 4
    
    return distance_dict

def compute_host_microbiome_similarities(ds,
                                         config,
                                         treatment,
                                         within_studies=False):
    all_hosts = list(ds.metadata_df["HostSpecific"].unique())
    tax_df = ds.taxonomy_counts_df
    similarity_dict = defaultdict(dict)

    for host1 in all_hosts:
        for host2 in all_hosts:
            if host1 == host2:
                continue

            mdf = ds.metadata_df

            mask1 = (mdf["HostSpecific"] == host1) & (mdf["Treatment"] == treatment)
            mask2 = (mdf["HostSpecific"] == host2) & (mdf["Treatment"] == treatment)

            if not within_studies:
                common_studies = []
                for study in config["drought_studies"]:
                    hosts_map = config[study]['hosts_ungrouped']
                    if host1 in hosts_map.values() and host2 in hosts_map.values():
                        common_studies.append(study)

                if common_studies:
                    excl_mask = mdf["StudyID"].isin(common_studies) & (
                        mdf["HostSpecific"].isin([host1, host2]) &
                        (mdf["Treatment"] == treatment)
                    )
                    mask1 = mask1 & ~excl_mask
                    mask2 = mask2 & ~excl_mask

            if not mask1.any() or not mask2.any():
                similarity_dict[host1][host2] = np.nan
                continue

            counts1 = tax_df.loc[mask1, :]
            counts2 = tax_df.loc[mask2, :]

            pres1 = (counts1 > 0).any(axis=0)
            pres2 = (counts2 > 0).any(axis=0)

            intersection = np.logical_and(pres1, pres2).sum()
            similarity_dict[host1][host2] = intersection

    return similarity_dict
