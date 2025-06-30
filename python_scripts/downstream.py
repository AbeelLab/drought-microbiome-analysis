import argparse
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle
from scipy.stats import spearmanr
import yaml

from analysis import run_limma_diff_abundance
from dataset import Dataset
from plotting import make_volcano_plot, plot_correlations, draw_network, plot_lfc_diff_abundance
from taxonomy_utils import get_last_taxonomic_level
from utils import biom_to_tsv, process_metadata

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
    # Process metadata
    metadata_df = process_metadata(study_path,
                                   study,
                                   config[study])
            
    # Initialize Dataset
    dataset_name = f"{study}_l{level}"
    ds = Dataset(taxonomy_counts_df=taxonomy_counts_df,
                 metadata_df=metadata_df,
                 dataset_name=dataset_name,
                 data_path=config["data_path"],
                 is_filtered=False)
    ds.filter_features()
    print("Samples: ", len(taxonomy_counts_df))
    # Save filtered Dataset as .tsv
    ds.save_dataset()
    
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
        run_preprocessing(config["inoculum_studies"], [6], config["processed_inoculum_datasets"])

    with open(config["processed_drought_datasets"], 'rb') as handle:
        processed_drought_datasets = pickle.load(handle)
    with open(config["processed_inoculum_datasets"], 'rb') as handle:
        processed_inoculum_datasets = pickle.load(handle)        
    
    for level, level_name in [[6, "genus"]]:
        # Merge drought datasets
        dataset_name = f"merged_l{level}"
        list_of_datasets = [processed_drought_datasets[level][study]
                            for study in config["drought_studies"]]
        merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                config["data_path"],
                                                merged_dataset_name=dataset_name)
        # Save merged dataset
        merged_dataset.save_dataset()

        # ---> Downstream analysis
        # Find drought signature
        correction = "holm"
        if not os.path.exists(f"drought_signature_{correction}.pkl"):
            drought_signature = run_limma_diff_abundance(merged_dataset,
                                                         p_val=0.001,
                                                         correction=correction)
            with open(f"drought_signature_{correction}.pkl", 'wb') as handle:
                pickle.dump(drought_signature,
                            handle,
                            protocol=pickle.HIGHEST_PROTOCOL)
        else:
            with open(f"drought_signature_{correction}.pkl", 'rb') as handle:
                drought_signature = pickle.load(handle)

        # Filter drought signature
        drought_signature = {taxon: {'logFC': drought_signature[taxon]['logFC'],
                                     'adj.P.Val': drought_signature[taxon]['adj.P.Val']}
                             for taxon in drought_signature
                             if abs(drought_signature[taxon]['logFC']) >= 0.25}
        print("Drought signature size:", len(drought_signature))

        plot_lfc_diff_abundance(drought_signature,
                                level,
                                config)

        # Further process inoculum studies
        for treatment in ["Control", "Drought"]:
            per_treatment_datasets = dict()
            for study, inoculum_dataset in processed_inoculum_datasets[level].items():
                ds = deepcopy(inoculum_dataset)
                # Remove rows that are not Plant-associated (i.e. inoculum soil samples)
                ds.filter_rows(lambda df: df['Sample type'] == 'Plant-associated')
                # Keep only Drought samples
                ds.filter_rows(lambda df: df['Treatment'] == treatment)
                # Remove sterilized controls
                # (only compare wet and drought-legacy inocula)
                ds.filter_rows(lambda df: df['InoculumSubtypes'] != "Control")

                print(f"[INFO] Re-filtering features after removing samples")
                # After removing samples, re-do feature filtering (less strict)
                ds.filter_features()

                if level not in per_treatment_datasets:
                    per_treatment_datasets[level] = {study: ds}
                else:
                    per_treatment_datasets[level][study] = ds
               
            inoculum_studies = {"positive": config["positive_inoculum_studies"],
                                "negative": config["negative_inoculum_studies"]}
            merged_inoculum_datasets_by_outcome = dict()
            for outcome, outcome_study_list in inoculum_studies.items():
                # Merge per outcome type
                dataset_name = f"merged_inoculum_{outcome}_l{level}"
                list_of_datasets = [per_treatment_datasets[level][study]
                                    for study in outcome_study_list]
                merged_dataset = Dataset.merge_datasets(list_of_datasets,
                                                        config["data_path"],
                                                        merged_dataset_name=dataset_name)
                merged_dataset.save_dataset()
                merged_inoculum_datasets_by_outcome[outcome] = merged_dataset

            # Perform differential abundance by inoculum type
            da_inoculum_by_outcome = dict()
            for outcome, dataset in merged_inoculum_datasets_by_outcome.items():
                # p_val set to 1 so we return all for plotting
                da_results = run_limma_diff_abundance(dataset,
                                                      p_val=1,
                                                      formula="Inoculum + Study",
                                                      variable="Inoculum",
                                                      base="Control",
                                                      condition="DroughtLegacy")
                da_inoculum_by_outcome[outcome] = da_results

            # Make volcano plot per outcome, highlight drought signature
            for outcome, da_results in da_inoculum_by_outcome.items():
                save_as = os.path.join(config["plotting_dir"],
                                       f"volcano_{outcome}_l{level}_{treatment.lower()}.svg")
                n = len(merged_inoculum_datasets_by_outcome[outcome].get_counts_features_columns())
                make_volcano_plot(da_inoculum_by_outcome[outcome],
                                  drought_signature,
                                  save_as,
                                  total_num_features=n)

            # Case study: outliers for Drought, positive outcome
            if treatment == "Drought":
                # Per study analysis
                x = da_inoculum_by_outcome["positive"]
                outlier_taxa = [taxa for taxa in x
                                if (abs(x[taxa]["logFC"]) >= 2 and x[taxa]["adj.P.Val"] <= 0.05)]
                outlier_taxa_signature = [taxa for taxa in outlier_taxa
                                          if taxa in drought_signature]
                print("Outlier drought signature taxa:", len(outlier_taxa_signature))

                # --- Check if these are confirmed by the two studies individually
                for study in inoculum_studies["positive"]:
                    print("Study: ", study)
                    ds = per_treatment_datasets[level][study]
                    da_per_study = run_limma_diff_abundance(ds,
                                                            p_val=0.05,
                                                            formula="Inoculum",
                                                            variable="Inoculum",
                                                            base="Control",
                                                            condition="DroughtLegacy")

                    common = set(da_per_study).intersection(set(outlier_taxa_signature))
                    print("Common: ", len(common))
                    print(common)

                # --- Correlation analysis for successful inoculation
                save_as = os.path.join(config["plotting_dir"],
                                       f"correlations_with_signature_positive_l{level}.svg")
                merged_inoculum_datasets_by_outcome["positive"].apply_clr()
                assert merged_inoculum_datasets_by_outcome["positive"].is_clr_transformed
                plot_correlations(da_inoculum_by_outcome["positive"],
                                  drought_signature,
                                  merged_inoculum_datasets_by_outcome["positive"],
                                  save_as)

                # --- Heatmap
                matrix = [[0, 0], [0, 0]]
                for taxon in da_inoculum_by_outcome["positive"]:
                    if taxon in drought_signature:
                        if drought_signature[taxon]['logFC'] > 0 and da_inoculum_by_outcome["positive"][taxon]['logFC'] > 0:
                            matrix[1][1] += 1
                        if drought_signature[taxon]['logFC'] > 0 and da_inoculum_by_outcome["positive"][taxon]['logFC'] < 0:
                            matrix[1][0] += 1
                        if drought_signature[taxon]['logFC'] < 0 and da_inoculum_by_outcome["positive"][taxon]['logFC'] > 0:
                            matrix[0][1] += 1
                        if drought_signature[taxon]['logFC'] < 0 and da_inoculum_by_outcome["positive"][taxon]['logFC'] < 0:
                            matrix[0][0] += 1
                print(matrix)
                
                depleted = 0
                enriched = 0
                for taxon in da_inoculum_by_outcome["positive"]:
                    if da_inoculum_by_outcome["positive"][taxon]['logFC'] > 0:
                        enriched +=1
                    else:
                        depleted += 1
                print("Enriched: ", enriched)
                print("Depleted: ", depleted)
                
                print(drought_signature)

                # --- Bar plot
                data = dict()
                for taxon in outlier_taxa_signature:
                    data[get_last_taxonomic_level(taxon)] = {'inoculum': da_inoculum_by_outcome["positive"][taxon]['logFC'],
                                                             'drought signature': drought_signature[taxon]['logFC']}

                # horizontal bar plot
                # per taxon: two bars
                # one for inoculum, one for drought signature
                taxa = list(data.keys())
                inoculum_vals = [data[t]['inoculum'] for t in taxa]
                drought_vals = [data[t]['drought signature'] for t in taxa]

                # Bar positions
                ind = np.arange(len(taxa))
                width = 0.35

                fig, ax = plt.subplots(figsize=(3, len(taxa)))
                ax.barh(ind - width/2, inoculum_vals, height=width, label='Inoculum', color='steelblue')
                ax.barh(ind + width/2, drought_vals, height=width, label='Drought signature', color='salmon')
                    
                ax.set_yticks(ind)
                ax.set_yticklabels(taxa)
                ax.axvline(0, color='grey', linewidth=0.8)
                ax.set_xlabel('log2 Fold Change')
                ax.set_title('Outlier Taxa – Positive Outcome')
                ax.legend()
                    
                plt.tight_layout()

                save_as = os.path.join(config["plotting_dir"],
                                           f"bar_outlier_positive.svg")
                plt.savefig(save_as)
                plt.close()
                    

                # --- Network analysis
                for inoculation_status in ["DroughtLegacy", "Control"]:                    
                    # Make edge list
                    edge_list = []
                    ds = deepcopy(merged_inoculum_datasets_by_outcome["positive"])
                    ds.filter_rows(lambda df: df['Inoculum'] == inoculation_status)
                    for taxa1 in outlier_taxa_signature:
                        for taxa2 in outlier_taxa_signature:
                            if taxa1 != taxa2:
                                res = spearmanr(ds.get_counts_feature(taxa1),
                                                ds.get_counts_feature(taxa2))
                                if abs(res.statistic) > 0.5:
                                    edge_list.append((get_last_taxonomic_level(taxa1),
                                                      get_last_taxonomic_level(taxa2),
                                                      res.statistic))
                                else:
                                    edge_list.append((get_last_taxonomic_level(taxa1),
                                                      get_last_taxonomic_level(taxa2),
                                                      0))
                                        
                    edge_df = pd.DataFrame(edge_list,
                                           columns=['source', 'target', 'weight'])
                
                    node_colors = dict()
                    for node in outlier_taxa:
                        if node in drought_signature and drought_signature[node]["logFC"] < 0:
                            node_colors[get_last_taxonomic_level(node)] = "black"
                        elif node in drought_signature and drought_signature[node]["logFC"] >= 0:
                            node_colors[get_last_taxonomic_level(node)] = "#fc8135ff"
                        else:
                            node_colors[get_last_taxonomic_level(node)] = "gray"

                    save_as = os.path.join(config["plotting_dir"],
                                       f"network_positive_l{level}_{inoculation_status.lower()}.svg")
                    draw_network(edge_df,
                                 node_colors,
                                 save_as)

                for taxa in outlier_taxa_signature:
                    print(treatment)
                    print(taxa)
                    print(x[taxa])
                    print(drought_signature[taxa])

if __name__ == "__main__":
    main()
