import os
import pickle
import yaml

from collections import OrderedDict
from copy import deepcopy
from dataset import Dataset
from plotting import plot_pcoa
from process_per_study import common_processing
from utils import biom_to_tsv, log_statistics, log_features_per_batch, get_last_taxonomic_level

with open('config.yml') as f:
    config = yaml.safe_load(f)

def create_study_dataset(study, level):
    print(f"[INFO] Processing study: {study}...")
    study_path = os.path.join(config["data_path"],
                              study)

    # Convert .biom files produced by QIIME
    comp_dir = os.path.join(study_path,
                            f"composition_table_l{level}")
    taxonomy_counts_df = biom_to_tsv(comp_dir)
    kept_samples = taxonomy_counts_df.index

    # Process metadata
    custom_processing = None
    if "custom_processing" in config[study]:
        custom_processing = config[study]["custom_processing"]
    metadata_df = common_processing(study_path,
                                    study,
                                    config[study],
                                    custom_processing,
                                    kept_samples)
        
    # Initialize Dataset
    level_name = config["levels"][level]
    dataset_name = f"{level_name}"
    ds = Dataset(taxonomy_counts_df=taxonomy_counts_df,
                 metadata_df=metadata_df,
                 dataset_name=dataset_name,
                 data_path=os.path.join(config["processed_data_path"], study),
                 study=study)

    ds.remove_features_not_assigned_at_last_level()
    
    # Shouldn't happen, but as a sanity check
    ds.remove_zero_features()
    if level == 6:
        ds.save_dataset(save_metadata=True)
    else:
        ds.save_dataset()

    # Log sample statistics
    if level == 6:
        ds.log_sample_statistics()
        ds.log_sparsity_and_feature_statistics()
    
    return ds

def run_preprocessing(studies, save_as):
    processed_datasets = dict()

    for level in config["levels"]:
        print(f"---\n[INFO] Level: {level}")
        processed_datasets[level] = dict()
        for study in studies:
            ds = create_study_dataset(study, level)
            processed_datasets[level][study] = ds
    
    # Pickle Dataset dictionary
    with open(save_as, 'wb') as handle:
        pickle.dump(processed_datasets,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL)

def main():
    #run_preprocessing(config["all_studies"], config["processed_datasets"])

    with open(config["processed_datasets"], 'rb') as handle:
        processed_datasets = pickle.load(handle)

    for level in [2]:#config["levels"]:
        for compartment in ["all", "Rhizosphere", "Endosphere", "Bulk soil"]:
            list_of_datasets = []
            drought_datasets = []
            for study in config["all_studies"]:
                ds = deepcopy(processed_datasets[level][study])
                if compartment != "all":
                    ds.filter_samples_based_on_condition(lambda df: df["RootCompartment"] == compartment)
                    # Remove features not occuring in this compartment
                    ds.remove_zero_features()
                if len(ds.taxonomy_counts_df):
                    list_of_datasets.append(ds)
                    if study in config["drought_studies"]:
                        drought_datasets.append(ds)

            level_name = config["levels"][level]
            merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                    os.path.join(config["processed_data_path"],
                                                                 compartment.lower().replace(' ', '_')),
                                                    merged_dataset_name=f"{level_name}",
                                                    method="union")

            # For all compartments, no need for a sub-folder
            if compartment == "all":
                merged_dataset.data_path = config["processed_data_path"]
                if level == 6:
                    #Save all metadata only once
                    merged_dataset.save_dataset(save_metadata=True)

                if level == 2:
                    # Batch-correct phylum level counts across all compartments
                    merged_dataset_copy = deepcopy(merged_dataset)
                    merged_dataset_copy.apply_mmuphin_be_correction(covariates=["Treatment", "RootCompartment"])
                    merged_dataset_copy.remove_zero_features()
                    merged_dataset_copy.save_dataset()
                else:
                    merged_dataset.save_dataset()
            else:
                merged_dataset.save_dataset()

            if compartment == "all":
                print("[INFO] # hosts: ", len(merged_dataset.metadata_df["HostSpecific"].unique()))
                with open("../hosts.txt", "w") as f:
                    for host in merged_dataset.metadata_df["HostSpecific"].unique():
                        f.write(str(host) + "\n")

            print("[INFO] Merged features: ",
                  len(merged_dataset.taxonomy_counts_df.columns))
                
                
            # For genus-level and drought studies,
            # additionally merge only features that appear in 30% of datasets
            if level == 6:
                if compartment == "all":
                    data_path = config["processed_data_path"]
                else:
                    data_path = os.path.join(config["processed_data_path"],
                                             compartment.lower().replace(' ', '_'))
                merged_dataset_drought = Dataset.merge_datasets(drought_datasets,
                                                                data_path,
                                                                merged_dataset_name=f"{level_name}_for_signature",
                                                                method="intersection")

                # Remove uncultured taxa
                # So we set min_prevalence to 0
                merged_dataset_drought.filter_features(min_prevalence=0)
                merged_dataset_drought.save_dataset()

                merged_dataset_drought.apply_mmuphin_be_correction()
                assert "batch_corrected" in merged_dataset_drought.dataset_name
                merged_dataset_drought.remove_zero_features()
                merged_dataset_drought.save_dataset()

if __name__ == "__main__":
    main()
