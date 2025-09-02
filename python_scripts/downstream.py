import argparse
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
import yaml

from collections import OrderedDict
from dataset import Dataset
from plotting import make_volcano_plot, plot_correlations, draw_network, plot_lfc_diff_abundance, plot_pcoa, plot_heatmap, plot_corr, plot_taxonomy, plot_num_samples
from taxonomy_utils import get_last_taxonomic_level
from utils import biom_to_tsv, log_statistics, log_features_per_batch, compute_host_phylogenetic_distances, compute_host_microbiome_similarities, get_core_microbiome
from process_per_study import common_processing

with open('config.yml') as f:
    config = yaml.safe_load(f)

def create_and_filter_study_dataset(study, level):
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

    # Log #features and filter
    log_statistics(study,
                   "Initial #features that are not Chloroplast;Mitochondria;Eukaryota;Unassigned;Unclassified",
                   int(len(ds.taxonomy_counts_df.columns)),
                   "../data/features_log.pkl")
    # ds.filter_features()
    # log_statistics(study,
    #                "After abundance and prevalence filtering",
    #                int(len(ds.taxonomy_counts_df.columns)),
    #                "../data/features_log.pkl")

    log_statistics(study,
                   "Sparsity",
                   (ds.taxonomy_counts_df == 0).values.sum() / ds.taxonomy_counts_df.size * 100,
                   "../data/features_log.pkl")
    log_statistics(study,
                   "Average summed sample abundance",
                   ds.taxonomy_counts_df.sum(axis=1).mean(),
                   "../data/features_log.pkl")
    
    print("Samples: ", len(ds.taxonomy_counts_df))
    # Save filtered Dataset as .tsv
    ds.save_dataset()

    # Log #samples
    log_statistics(study,
                  "After filtering based on #ASVs and #taxa and removing hosts with too few samples",
                  len(ds.taxonomy_counts_df))
    
    return ds

def run_preprocessing(studies, levels, save_as):
    processed_datasets = dict()

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
    parser = argparse.ArgumentParser()
    parser.add_argument("--preprocess",
                        action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    preprocess = args.preprocess

    if preprocess:
        run_preprocessing(config["drought_studies"], [6], config["processed_drought_datasets"])
        # run_preprocessing(config["inoculum_studies"], [6], config["processed_inoculum_datasets"])

    with open(config["processed_drought_datasets"], 'rb') as handle:
        processed_drought_datasets = pickle.load(handle)
    # with open(config["processed_inoculum_datasets"], 'rb') as handle:
    #     processed_inoculum_datasets = pickle.load(handle)        

    
    for level, level_name in [[6, "genus"]]:
        # Merge everything
        dataset_name = f"merged_l{level}"
        list_of_datasets = [processed_drought_datasets[level][study]
                            for study in config["drought_studies"]]
        merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                config["data_path"],
                                                merged_dataset_name=dataset_name)
        merged_dataset.save_dataset()

        host_to_studies = OrderedDict()
        for ds in list_of_datasets:
            study_name = ds.metadata_df["StudyName"][0] 
            for host in ds.metadata_df["Host"]:
                host_to_studies.setdefault(host, set()).add(study_name)
        for host, studies in host_to_studies.items(): host_to_studies[host] = sorted(studies)

        ordered_ticks = []
        for host in host_to_studies:
            studies = list(host_to_studies[host])
            studies.sort()
            for study in studies:
                ordered_ticks.append([host, study])

        print(ordered_ticks)

        # Merge drought datasets per compartment
        for compartment in ["Rhizosphere", "Endosphere", "Bulk soil"]:
            print(compartment)
            dataset_name = f"merged_l{level}_{compartment.lower().replace(' ', '_')}"

            list_of_datasets = []
            for study in config["drought_studies"]:
                ds = deepcopy(processed_drought_datasets[level][study])
                ds.filter_rows(lambda df: df["RootCompartment"] == compartment)
                ds.dataset_name += "_" + compartment

                if len(ds.taxonomy_counts_df): list_of_datasets.append(ds)
                    
            merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                    config["data_path"],
                                                    merged_dataset_name=dataset_name)
            
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
            
            # Plot dataset statistics
            core = get_core_microbiome(merged_dataset)
            if compartment != "Bulk soil":
                plot_taxonomy(merged_dataset,
                              [f"{x} | {y}" for x,y in ordered_ticks
                               if x != "Soil"],
                              config["colors"]["Phyla"],
                              f"../data/plots/taxonomy_{compartment}")
                plot_num_samples(merged_dataset,
                                 [f"{x} | {y}" for x,y in ordered_ticks
                               if x != "Soil"],
                                 f"../data/plots/samples_{compartment}")
            else:
                plot_taxonomy(merged_dataset,
                              [f"{x} | {y}" for x,y in ordered_ticks
                               if x == "Soil"],
                              config["colors"]["Phyla"],
                              f"../data/plots/taxonomy_{compartment}")
                plot_num_samples(merged_dataset,
                                 [f"{x} | {y}" for x,y in ordered_ticks
                                  if x == "Soil"],
                                 f"../data/plots/samples_{compartment}")

            if compartment != "Bulk soil":
                for within_studies in [False, True]:
                    for treatment, alpha in [("Drought", 0.9), ("Control", 0.6)]:
                        phy_distances = compute_host_phylogenetic_distances()
                        similarities = compute_host_microbiome_similarities(merged_dataset,
                                                                            config,
                                                                            treatment,
                                                                            within_studies)
                
                        suffix = "_within_studies" if within_studies else ""
                        # Plot each pair as a dot (skip NaNs)
                        # Save as both svg and png
                        # x-axis: microbiome similarity
                        # y-axis: phylogenetic distance
                        plot_corr(phy_distances,
                                  similarities,
                                  dot_color=config["colors"]["RootCompartment"][compartment],
                                  alpha=alpha,
                                  save_as=f"../data/plots/corr_{compartment}_{treatment}{suffix}")

if __name__ == "__main__":
    main()
