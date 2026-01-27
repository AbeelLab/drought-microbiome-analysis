import math
import numpy as np
import pandas as pd
import pickle as pkl
import os

from collections import Counter
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from taxonomy_utils import get_last_taxonomic_level
from utils import log_statistics, get_last_taxonomic_level

class Dataset:
    def __init__(self,
                 taxonomy_counts_df,
                 metadata_df,
                 dataset_name,
                 data_path,
                 study=None):
        # Keep only common indices from both dataframes
        # (we want to keep samples with both metadata and taxonomic profiles)
        self.taxonomy_counts_df = taxonomy_counts_df
        self.metadata_df = metadata_df
        self.update_based_on_sample_intersection()

        assert len(self.taxonomy_counts_df) == len(self.metadata_df)

        # For saving
        self.dataset_name = dataset_name
        self.data_path = data_path
        self.study = study
        
        self.is_clr_transformed = False
        self.is_batch_corrected = False

    def update_based_on_sample_intersection(self):
        left = self.taxonomy_counts_df.index
        right = self.metadata_df.index
        # Do not sort: keep order from left dataframe
        common_samples = left.intersection(right,
                                           sort=False)
        self.taxonomy_counts_df = self.taxonomy_counts_df.loc[common_samples]
        self.metadata_df = self.metadata_df.loc[common_samples]

    def get_metadata(self):
        return self.metadata_df

    def get_metadata_columns(self):
        return self.metadata_df.columns.tolist()

    def get_counts_features(self):
        return self.taxonomy_counts_df

    def get_counts_features_columns(self):
        return self.taxonomy_counts_df.columns.tolist()

    def apply_clr(self,
                  pseudocount=1e-6):
        self.taxonomy_counts_df = self.taxonomy_counts_df + pseudocount
        log_df = np.log(self.taxonomy_counts_df)
        clr_df = log_df.subtract(log_df.mean(axis=1), axis=0)
        self.taxonomy_counts_df = clr_df
        self.is_clr_transformed = True
        self.dataset_name += "_clr"

    # Based on the MMuPHin user tutorial: https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html
    def apply_mmuphin_be_correction(self, covariates=["Treatment"]):
        pandas2ri.activate()
        mmuphin = importr('MMUPHin')

        counts_df = self.get_counts_features()
        metadata_df = self.get_metadata()
        assert not counts_df.isna().values.any(), "Count table includes NaNs"

        counts_r = pandas2ri.py2rpy(counts_df.T)
        metadata_r = pandas2ri.py2rpy(metadata_df)

        robjects.globalenv['counts'] = counts_r
        robjects.globalenv['metadata'] = metadata_r

        if len(covariates) == 0:
            covariate_arg = "NULL"
        elif len(covariates) == 1:
            covariate_arg = f'"{covariates[0]}"'
        else:
            covariate_arg = "c(" + ", ".join([f'"{c}"' for c in covariates]) + ")"

        robjects.r(f'''
        library(MMUPHin)
        fit_adjust_batch <- adjust_batch(
        feature_abd = counts,
        batch = "StudyID",
        covariates = {covariate_arg},
        data = metadata,
        control = list(verbose = TRUE)
        )

        adj_data <- as.data.frame(fit_adjust_batch$feature_abd_adj)
        ''')

        counts_batch_corrected = robjects.globalenv['adj_data']
        counts_batch_corrected = pandas2ri.rpy2py(counts_batch_corrected).T
        
        self.taxonomy_counts_df = counts_batch_corrected
        self.dataset_name += "_batch_corrected"
        self.is_batch_corrected = True

        self.update_based_on_sample_intersection()
        
        
    # Normalize samples such that the features sum up to total_sum
    def get_normalized_features(self, total_sum=100):
        df = self.taxonomy_counts_df.copy()
        taxonomy_normalized_df = df.div(df.sum(axis=1), axis=0) * total_sum
        return taxonomy_normalized_df

    def get_prevalence(self, feature):
        df = self.taxonomy_counts_df
        non_zero_count = (df[feature] > 0).sum()
        total_samples = df.shape[0]
    
        prevalence_percent = (non_zero_count / total_samples) * 100
        return prevalence_percent
         
    # Default: appears in at least 5% of samples
    # This function also removes uncultured taxa 
    def filter_features(self,
                        min_prevalence=0.05,
                        min_abundance=0,
                        total_sum=100):
        all_features = self.get_counts_features_columns()
        
        # Remove uncultured taxa
        to_keep = [feature for feature in all_features
                   if "uncultured" not in feature
                   and "Uncultured" not in feature]
        self.taxonomy_counts_df = self.taxonomy_counts_df[to_keep]

        # Normalize to enable abundance/prevalence filtering
        taxonomy_normalized_df = self.get_normalized_features(total_sum)

        # Scale abundance based on total_sum
        min_abundance *= total_sum
        
        min_sample_count = int(min_prevalence * len(taxonomy_normalized_df))
        to_keep = [col for col in to_keep
                   if (taxonomy_normalized_df[col] > min_abundance).sum() >= min_sample_count]

        # Update taxonomy dataframe
        self.taxonomy_counts_df = self.taxonomy_counts_df[to_keep]
        

    # Filter rows of the metadata based on a condition
    # Condition should specify which rows to keep
    def filter_samples_based_on_condition(self, condition):        
        self.metadata_df = self.metadata_df[condition(self.metadata_df)]
        self.update_based_on_sample_intersection()

    def save_dataset(self,
                     save_metadata=False,
                     metadata_file="metadata.tsv"):
        # Metadata doesn't have to be per taxonomic level, so don't add the dataset name
        if save_metadata:
            metadata_file = os.path.join(self.data_path, metadata_file)
            self.metadata_df.fillna("").to_csv(metadata_file,
                                               sep='\t',
                                               index=True,
                                               index_label="#SampleID")

        taxonomy_counts_file = "counts_" + self.dataset_name + ".tsv"
        taxonomy_counts_file = os.path.join(self.data_path, taxonomy_counts_file)
        self.taxonomy_counts_df.to_csv(taxonomy_counts_file,
                                       sep = '\t',
                                       index = True,
                                       index_label = "#SampleID")
        # Save transpose for .biom convesion
        transposed_file = "counts_" + self.dataset_name + ".transposed.tsv"
        transposed_file = os.path.join(self.data_path, transposed_file)
        self.taxonomy_counts_df.T.to_csv(transposed_file,
                                         sep = '\t',
                                         index = True,
                                         index_label = "#SampleID")

        with open(os.path.join(self.data_path, self.dataset_name + ".pkl"), 'wb') as handle:
            pkl.dump(self,
                     handle,
                     protocol=pkl.HIGHEST_PROTOCOL)

        return metadata_file, taxonomy_counts_file

    @staticmethod
    # Return a dataset with the merged metadata and taxonomic features
    def merge_datasets(list_of_datasets,
                       data_path,
                       method="union",
                       merged_dataset_name="merged",
                       intersection_threshold=0.3):
        merged_metadata = pd.concat([ds.get_metadata()
                                     for ds in list_of_datasets],
                                    axis=0,
                                    sort=False)

        if method == "union":
            # For features that are not common, pad with 0's
            merged_taxonomy_counts = pd.concat([ds.get_counts_features()
                                                for ds in list_of_datasets],
                                               axis=0,
                                               sort=False).fillna(0)
        elif method == "intersection":
            n_datasets = len(list_of_datasets)

            if intersection_threshold < 1:
                min_datasets = math.ceil(intersection_threshold * n_datasets)
            else:
                min_datasets = intersection_threshold

            counter = Counter()
            for ds in list_of_datasets:
                counter.update(ds.get_counts_features().columns)

            features_kept = {feat for feat, cnt in counter.items() if cnt >= min_datasets}
            ordered_cols = sorted(features_kept)

            dfs = [ds.get_counts_features().reindex(columns=ordered_cols, fill_value=0)
                   for ds in list_of_datasets]
            merged_taxonomy_counts = pd.concat(dfs, axis=0, sort=False)

        merged = Dataset(taxonomy_counts_df=merged_taxonomy_counts,
                         metadata_df=merged_metadata,
                         dataset_name=merged_dataset_name,
                         data_path=data_path)

        assert len(merged_metadata) == len(merged_taxonomy_counts)

        return merged

    def get_core(self, min_prevalence):
        all_features = self.get_counts_features_columns()
        n_samples = len(self.taxonomy_counts_df)

        min_sample_count = int(min_prevalence * n_samples)

        prevalence = {col: (self.taxonomy_counts_df[col] > 0).sum() / n_samples
                      for col in all_features}

        core_features = {taxon for taxon, prev in prevalence.items()
                         if prev >= min_prevalence}

        return core_features, prevalence

    
    def remove_zero_features(self):
        all_features = self.get_counts_features_columns()
        nonzero_features = self.taxonomy_counts_df.columns[self.taxonomy_counts_df.sum(axis=0) > 0].tolist()
        self.taxonomy_counts_df = self.taxonomy_counts_df[nonzero_features]

        
    def remove_features_not_assigned_at_last_level(self):
        all_features = self.get_counts_features_columns()
        to_keep = [feature for feature in all_features
                   if get_last_taxonomic_level(feature) is not None]
        self.taxonomy_counts_df = self.taxonomy_counts_df[to_keep]
        

    def log_sparsity_and_feature_statistics(self):
        df = self.taxonomy_counts_df
        name = self.study

        n_samples, n_features = df.shape
        total_elements = n_samples * n_features if (n_samples > 0 and n_features > 0) else 0

        zero_count = int((df == 0).sum().sum()) if total_elements > 0 else 0
        sparsity = zero_count / total_elements * 100 if total_elements > 0 else float('nan')

        # Min non-zero count
        stacked_nonzero = df.stack()[df.stack() != 0] if total_elements > 0 else pd.Series(dtype=float)
        min_nonzero = float(stacked_nonzero.min())

        # Max count
        max_value = float(df.max().max())

        # Mean count per sample
        mean_count_per_sample = float(df.sum(axis=1).mean()) if n_samples > 0 else float('nan')

        # Mean count per feature
        mean_count_per_feature = float(df.sum(axis=0).mean()) if n_features > 0 else float('nan')
        
        # Log metrics
        logfile = "../features_log.pkl"
        log_statistics(name, "sparsity", float(sparsity), logfile)
        log_statistics(name, "min_nonzero_count", min_nonzero, logfile)
        log_statistics(name, "max_count_value", max_value, logfile)
        log_statistics(name, "mean_total_count_per_sample", mean_count_per_sample, logfile)
        log_statistics(name, "mean_count_per_feature", mean_count_per_feature, logfile)
        log_statistics(name, "# features", n_features, logfile)
        log_statistics(name, "# samples", n_samples, logfile)

        # Print summary
        print(f"[INFO] {name}: samples={n_samples}, features={n_features}, sparsity={sparsity:.4f}, "
              f"min_nonzero={min_nonzero}, max={max_value}, mean_per_sample={mean_count_per_sample:.3f}")

    def log_sample_statistics(self):
        metadata = self.get_metadata()
        logfile = "../samples_log.pkl"
        name = self.study

        total_samples = len(metadata)

        treatment_drought = treatment_control = None
        inoc_drought = inoc_ref = inoc_sterile = None

        if "Treatment" in metadata.columns:
            treatment_counts = metadata["Treatment"].value_counts(dropna=False).to_dict()
            treatment_drought = treatment_counts.get("Drought", 0)
            treatment_control = treatment_counts.get("Control", 0)

            treatment_sum = treatment_drought + treatment_control
            if treatment_sum > 0:
                assert treatment_sum == total_samples, (
                    f"[ERROR] Treatment counts ({treatment_sum}) "
                    f"do not match total samples ({total_samples}) for {name}")
        else:
            print(f"[WARN] 'Treatment' column not found in metadata for {name}")

        if "Inoculum" in metadata.columns:
            inoc_counts = metadata["Inoculum"].value_counts(dropna=False).to_dict()
            inoc_drought = inoc_counts.get("Dry", 0)
            inoc_ref = inoc_counts.get("Reference", 0)
            inoc_sterile = inoc_counts.get("Sterile", 0)

            inoc_sum = inoc_drought + inoc_ref + inoc_sterile
            if inoc_sum > 0:
                assert inoc_sum == total_samples, (
                    f"[ERROR] Inoculum counts ({inoc_sum}) do not match total samples ({total_samples}) for {name}"
                )
        else:
            print(f"[WARN] 'Inoculum' column not found in metadata for {name}")

        log_statistics(name, "total_samples", total_samples, logfile)
        if treatment_drought is not None:
            log_statistics(name, "Treatment_Drought", treatment_drought, logfile)
            log_statistics(name, "Treatment_Control", treatment_control, logfile)
        if inoc_drought is not None:
            log_statistics(name, "Inoculum_Dry", inoc_drought, logfile)
            log_statistics(name, "Inoculum_Reference", inoc_ref, logfile)
            log_statistics(name, "Inoculum_Sterile", inoc_sterile, logfile)

        # Print summary 
        print(f"[INFO] {name}: total_samples={total_samples}, "
              f"Treatment(Drought={treatment_drought}, Control={treatment_control}), "
              f"Inoculum(Dry={inoc_drought}, Reference={inoc_ref}, Sterile={inoc_sterile})")
