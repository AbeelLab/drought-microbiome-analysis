import os
import pickle
import yaml

from collections import OrderedDict
from copy import deepcopy
from dataset import Dataset
from plotting import plot_pcoa
from process_per_study import common_processing
from utils import biom_to_tsv, log_statistics, log_features_per_batch

with open('config.yml') as f:
    config = yaml.safe_load(f)

def create_and_filter_study_dataset(study,
                                    level):
    print(f"[INFO] Study: {study}")
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
        custom_processing = config[study]['custom_processing']
    metadata_df = common_processing(study_path,
                                    study,
                                    config[study],
                                    custom_processing,
                                    kept_samples)
            
    # Initialize Dataset
    dataset_name = f"{study}_l{level}"
    ds = Dataset(taxonomy_counts_df=taxonomy_counts_df,
                 metadata_df=metadata_df,
                 dataset_name=dataset_name,
                 data_path=config["data_path"],
                 is_filtered=False)

    ds.save_dataset()

    # Log #features and filter unculured, archaea
    log_statistics(study,
                   "Initial #features that are not Chloroplast;Mitochondria;Eukaryota;Unassigned;Unclassified",
                   int(len(ds.taxonomy_counts_df.columns)),
                   "../data/features_log.pkl")
    ds.filter_features(min_prevalence=0)
    log_statistics(study,
                   "After abundance and prevalence filtering",
                   int(len(ds.taxonomy_counts_df.columns)),
                   "../data/features_log.pkl")

    ds.filter_samples()
    
    ds.remove_zero_features()
    
    log_statistics(study,
                   "Sparsity",
                   (ds.taxonomy_counts_df == 0).values.sum() / ds.taxonomy_counts_df.size * 100,
                   "../data/features_log.pkl")
    log_statistics(study,
                   "Average summed sample abundance",
                   ds.taxonomy_counts_df.sum(axis=1).mean(),
                   "../data/features_log.pkl")
    log_statistics(study,
                   "After filtering based on absolute abundance",
                   len(ds.taxonomy_counts_df),
                   "../data/samples_log.pkl")
    
    print("Samples: ", len(ds.taxonomy_counts_df))
    # Save filtered Dataset as .tsv
    ds.save_dataset()

    # Log #samples
    log_statistics(study,
                  "After filtering based on #ASVs and #taxa and removing hosts with too few samples",
                  len(ds.taxonomy_counts_df))
    
    return ds

def run_preprocessing(studies, levels, save_as):
    processed_datasets = {"filtered": dict(),
                          "unfiltered": dict()}

    for level in levels:
        print(f"[INFO] Processing level: {level}")
        processed_datasets[level] = dict()
        
        for study in studies:
            ds = create_and_filter_study_dataset(study, level)
            processed_datasets[level][study] = ds
    
    # Pickle Dataset dictionary
    with open(save_as, 'wb') as handle:
        pickle.dump(processed_datasets,
                    handle,
                    protocol=pickle.HIGHEST_PROTOCOL)

def main():
    run_preprocessing(config["drought_studies"], [6], config["processed_drought_datasets"])
    run_preprocessing(config["inoculum_studies"], [6], config["processed_inoculum_datasets"])

    with open(config["processed_drought_datasets"], 'rb') as handle:
        processed_drought_datasets = pickle.load(handle)
    with open(config["processed_inoculum_datasets"], 'rb') as handle:
        processed_inoculum_datasets = pickle.load(handle)  
        
    
    for level, level_name in [[6, "genus"]]:
        dataset_name = f"merged_l{level}"
        list_of_datasets = [processed_drought_datasets[level][study]
                            for study in config["drought_studies"]]
        merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                    config["data_path"],
                                                    merged_dataset_name=dataset_name)
        merged_dataset.save_dataset()

        # Merge drought datasets per compartment
        for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
            dataset_name = f"merged_l{level}_{compartment.lower().replace(' ', '_')}"

            list_of_datasets = []
            for study in config["drought_studies"]:
                ds = deepcopy(processed_drought_datasets[level][study])
                ds.filter_rows(lambda df: df["RootCompartment"] == compartment)
                ds.dataset_name += "_" + compartment

                if len(ds.taxonomy_counts_df):
                    # Filter low-prevalence features
                    ds.filter_features(min_prevalence=0.05)
                    list_of_datasets.append(ds)

            for method in ["union", "intersection"]:
                suffix = ""
                if method == "intersection": suffix = "_intersection"
                merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                        config["data_path"],
                                                        merged_dataset_name=f"{dataset_name}{suffix}",
                                                        method=method)
                merged_dataset.remove_zero_features()
                    
                # Save merged dataset
                merged_dataset.save_dataset()
                plot_pcoa(merged_dataset,
                              os.path.join(config["plotting_dir"],
                                           merged_dataset.dataset_name))

                # Batch-correct merged dataset
                merged_dataset.apply_mmuphin_be_correction()
                merged_dataset.save_dataset()
                plot_pcoa(merged_dataset,
                          os.path.join(config["plotting_dir"],
                                       merged_dataset.dataset_name))
                log_statistics(merged_dataset.dataset_name,
                               "Union of features after BE correction filtering",
                               len(merged_dataset.taxonomy_counts_df.columns),
                               "../data/features_log.pkl")

            # Per study, log #features after batch effect removal (mmuphin filtering)
            log_features_per_batch(merged_dataset,
                                   compartment)

if __name__ == "__main__":
    main()
