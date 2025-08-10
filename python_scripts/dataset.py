import numpy as np
import pandas as pd
import pickle as pkl
import os

from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import Formula
from taxonomy_utils import get_last_taxonomic_level
from utils import log_statistics

class Dataset:
    def __init__(self,
                 taxonomy_counts_df,
                 metadata_df,
                 dataset_name,
                 data_path,
                 is_filtered):
        # Keep only common indices from both dataframes
        # (we want to keep samples with both metadata and taxonomic profiles)
        self.taxonomy_counts_df = taxonomy_counts_df
        self.metadata_df = metadata_df
        self.update_based_on_sample_intersection()

        assert len(self.taxonomy_counts_df) == len(self.metadata_df)

        # For saving
        self.dataset_name = dataset_name
        self.data_path = data_path
        
        self.is_filtered = is_filtered
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

    def get_counts_feature(self, feature):
        return self.taxonomy_counts_df[feature]

    def apply_clr(self,
                  pseudocount=1e-6):
        self.taxonomy_counts_df = self.taxonomy_counts_df + pseudocount
        log_df = np.log(self.taxonomy_counts_df)
        clr_df = log_df.subtract(log_df.mean(axis=1), axis=0)
        self.taxonomy_counts_df = clr_df
        self.is_clr_transformed = True
        self.dataset_name += "_clr"

    # Based on the MMuPHin user tutorial: https://bioconductor.org/packages/release/bioc/vignettes/MMUPHin/inst/doc/MMUPHin.html
    def apply_mmuphin_be_correction(self):
        pandas2ri.activate()
        
        mmuphin = importr('MMUPHin')
        
        # According to documentation, MMUPHin can handle counts 
        counts_df = self.get_counts_features()
        metadata_df = self.get_metadata()
        assert not counts_df.isna().values.any(), "Count table includes NaNs"

        assert "merged" in self.dataset_name, "Cannot run BE correction on single-study dataset"
        
        # Samples are columns
        counts_r = pandas2ri.py2rpy(counts_df.T)
        # Samples are rows
        metadata_r = pandas2ri.py2rpy(metadata_df)
        robjects.globalenv['counts'] = counts_r
        robjects.globalenv['metadata'] = metadata_r
        
        robjects.r(f'''
        library(MMUPHin)
        fit_adjust_batch <- adjust_batch(feature_abd = counts,
                                         batch = "StudyID",
                                         covariates = "Treatment",
                                         data = metadata,
                                         control = list(verbose = TRUE))

        adj_data <- as.data.frame(fit_adjust_batch$feature_abd_adj)
        ''')

        counts_batch_corrected = robjects.globalenv['adj_data']
        counts_batch_corrected = pandas2ri.rpy2py(counts_batch_corrected).T
        self.taxonomy_counts_df = counts_batch_corrected
        self.dataset_name += "_batch_corrected"

        self.update_based_on_sample_intersection()

    def filter_features_based_on_batches(self,
                                         batch_column="StudyID",
                                         min_num_batches=2):
        assert "merged" in self.dataset_name, "Need a merged dataset for this filtering step"
        
        all_features = self.get_counts_features_columns()
        sample_to_batch = self.metadata_df[batch_column]
        feature_batches = {feature: set() for feature in all_features}

        for batch_id, sample_ids in sample_to_batch.groupby(sample_to_batch).groups.items():
            batch_df = self.taxonomy_counts_df.loc[sample_ids, all_features]
            present_features = batch_df.columns[(batch_df > 0).any(axis=0)]
            for feature in present_features:
                feature_batches[feature].add(batch_id)


            log_statistics(f"{batch_id}_{self.dataset_name}",
                           "After abundance and prevalence filtering",
                           int(len(present_features)),
                           "../data/features_log.pkl")
            

        to_keep = [feature for feature, batches in feature_batches.items()
                   if len(batches) >= min_num_batches]

        self.taxonomy_counts_df = self.taxonomy_counts_df[to_keep]

        is_filtered = True

        print(f"[INFO] After filtering to keep featrues that occur in at least {min_num_batches} studies")
        print(f"{len(to_keep)} features")

        log_statistics(self.dataset_name,
                       "Union of features after filtering",
                       len(to_keep),
                       "../data/features_log.pkl")
        
    # Normalize samples such that the features sum up to total_sum
    def get_normalized_features(self, total_sum=100):
        df = self.taxonomy_counts_df.copy()
        taxonomy_normalized_df = df.div(df.sum(axis=1), axis=0) * total_sum
        return taxonomy_normalized_df

    
    # Default: at least > 0.01% abundant in at least 5% of samples
    def filter_features(self,
                        min_prevalence=0.05,
                        min_abundance=0.00001,
                        total_sum=100):        
        # Normalize features
        # Remove taxa not assigned at the deepest level
        all_features = self.get_counts_features_columns()
        print(f"[INFO] Initial #features: {len(all_features)}")
        
        to_keep = [feature for feature in all_features
                   if get_last_taxonomic_level(feature) is not None]
        # Remove uncultured taxa
        to_keep = [feature for feature in to_keep
                   if "uncultured" not in feature
                   and "Uncultured" not in feature]
        self.taxonomy_counts_df = self.taxonomy_counts_df[to_keep]

        # Normalize to enable abundance/prevalence filtering
        taxonomy_normalized_df = self.get_normalized_features(total_sum)

        # Scale abundance based on total_sum
        min_abundance *= total_sum
        
        min_sample_count = int(min_prevalence * len(taxonomy_normalized_df))
        to_keep = [col for col in to_keep
                   if (taxonomy_normalized_df[col] >= min_abundance).sum() >= min_sample_count]

        print(f"[INFO] #features after filtering: {len(to_keep)}")

        # Update taxonomy dataframe
        self.taxonomy_counts_df = self.taxonomy_counts_df[to_keep]
        
        # Update filtering status
        is_filtered = True
        self.dataset_name += "_filtered"

        
    def save_dataset(self):
        metadata_file = "metadata_" + self.dataset_name + ".tsv"
        metadata_file = os.path.join(self.data_path, metadata_file)
        self.metadata_df.to_csv(metadata_file,
                                sep = '\t',
                                index = True)
        print(f"[INFO] Saved: {metadata_file}")

        taxonomy_counts_file = "counts_" + self.dataset_name + ".tsv"
        taxonomy_counts_file = os.path.join(self.data_path, taxonomy_counts_file)
        self.taxonomy_counts_df.to_csv(taxonomy_counts_file,
                                       sep = '\t',
                                       index = True)
        print(f"[INFO] Saved: {taxonomy_counts_file}")

        with open(self.dataset_name + ".pkl", 'wb') as handle:
            pkl.dump(processed_datasets,
                     handle,
                     protocol=pickle.HIGHEST_PROTOCOL)

        return metadata_file, taxonomy_counts_file

    # Filter rows of the metadata based on a condition
    # Condition should specify which rows to keep
    def filter_rows(self, condition):
        print(f"[INFO] Filtering samples for {self.dataset_name}")
        print(f"[INFO] Initial #samples: {len(self.metadata_df)}")
        
        self.metadata_df = self.metadata_df[condition(self.metadata_df)]
        self.update_based_on_sample_intersection()

        print(f"[INFO] #samples after filtering: {len(self.metadata_df)}")

    @staticmethod
    # Return a dataset with the merged metadata and taxonomic features
    def merge_datasets(list_of_datasets, data_path, merged_dataset_name="merged"):
        filtering_stats = {ds.is_filtered for ds in list_of_datasets}
        assert len(filtering_stats) == 1, (f"Should not merge datasets with mixed filtering states")

        is_filtered = filtering_stats.pop()

        merged_metadata = pd.concat([ds.get_metadata()
                                     for ds in list_of_datasets],
                                    axis=0,
                                    sort=False)
        # For features that are not common, pad with 0's
        merged_taxonomy_counts = pd.concat([ds.get_counts_features()
                                            for ds in list_of_datasets],
                                           axis=0,
                                           sort=False).fillna(0)

        merged = Dataset(taxonomy_counts_df=merged_taxonomy_counts,
                         metadata_df=merged_metadata,
                         dataset_name=merged_dataset_name,
                         data_path=data_path,
                         is_filtered=is_filtered)

        assert len(merged_metadata) == len(merged_taxonomy_counts)

        print(f"[INFO] Merged dataset {merged_dataset_name}")
        print(f"with {len(merged_taxonomy_counts.columns)} features")
        print(f"and {len(merged_metadata)} samples")

        log_statistics(merged.dataset_name,
                       "Union of features",
                       len(merged.taxonomy_counts_df.columns),
                       "../data/features_log.pkl")

        return merged

    def get_core(self,
                 min_prevalence=0.95):
        min_sample_count = int(min_prevalence * len(self.taxonomy_counts_df))
        all_features = self.get_counts_features_columns()

        # Get columns (features) that appear in at least  min_sample_count

        # Return features as set
