import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle as pkl
import yaml

from copy import deepcopy
from matplotlib_venn import venn2
from matplotlib.colors import LogNorm
import seaborn as sns
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.metrics import f1_score
from sklearn.model_selection import GridSearchCV
from utils import trim_taxonomy

with open('config.yml') as f:
    config = yaml.safe_load(f)

def train_and_evaluate_rf(drought_ds,
                          study_ds,
                          features):
    drought_X = drought_ds.taxonomy_counts_df[features]
    drought_y = drought_ds.metadata_df.loc[drought_X.index, "Treatment"]
    drought_study_ids = drought_ds.metadata_df.loc[drought_X.index, "StudyID"]

    study_X = study_ds.taxonomy_counts_df[features]
    study_y = study_ds.metadata_df.loc[study_X.index, "Treatment"]

    param_grid = {
        "n_estimators": [100, 250],
        "max_features": ["sqrt", "log2"],
        "min_samples_leaf": [1, 0.05, 0.1, 0.25]
    }

    base_clf = RandomForestClassifier(random_state=42, n_jobs=1)

    # Split folds by StudyID
    unique_studies = drought_study_ids.unique()
    cv_splits = []
    for study in unique_studies:
        val_mask = (drought_study_ids == study).to_numpy()
        train_mask = ~val_mask
        val_idx = np.where(val_mask)[0]
        train_idx = np.where(train_mask)[0]
        cv_splits.append((train_idx, val_idx))

    grid = GridSearchCV(estimator=base_clf,
                        param_grid=param_grid,
                        cv=cv_splits,
                        scoring="f1_weighted",
                        n_jobs=1,
                        refit=True)

    grid.fit(drought_X, drought_y)

    best_clf = grid.best_estimator_
    best_idx = grid.best_index_
    score_cv = grid.cv_results_["mean_test_score"][best_idx]
    cv_std = grid.cv_results_["std_test_score"][best_idx]

    print("CV score:", score_cv)

    y_pred_test = best_clf.predict(study_X)
    score_test = f1_score(study_y, y_pred_test, average="weighted")

    y_pred_train = best_clf.predict(drought_X)
    score_train = f1_score(drought_y, y_pred_train, average="weighted")

    print("Training score:", score_train)
    print("Test score:", score_test)

    dummy = DummyClassifier(strategy="most_frequent", random_state=42)
    dummy.fit(drought_X, drought_y)
    dummy_pred = dummy.predict(study_X)
    dummy_score = f1_score(study_y, dummy_pred, average="weighted")
    print("Dummy test score:", dummy_score)
    
    dummy_pred_train = dummy.predict(drought_X)
    dummy_score_train = f1_score(drought_y, dummy_pred_train, average="weighted")
    print("Dummy train score:", dummy_score_train) 

    return best_clf, score_test, score_train, dummy_score, dummy_score_train, score_cv, cv_std


def make_grouped_bar_plot(f1_scores,
                          dummies,
                          stds,
                          x_ticks,
                          save_as="../data/plots/bar_rfs"):
    x = np.arange(len(x_ticks))
    width = 0.1

    fig, ax = plt.subplots(figsize=(12, 6))

    # Train bars
    rects1 = ax.bar(x - width*2, f1_scores["all_train"], width,
                    color="white", edgecolor="black", label="All (train)")
    rects2 = ax.bar(x - width, f1_scores["signature_train"], width,
                    color="black", label="Signature (train)")
    ax.scatter(x - width*2, dummies["all_train"], color="red", zorder=3)
    ax.scatter(x - width, dummies["signature_train"], color="red", zorder=3)

    # CV bars (with stds)
    rects3 = ax.bar(x, f1_scores["all_cv"], width,
                    yerr=stds["all_cv"], capsize=5,
                    color="white", edgecolor="black", label="All (CV)")
    rects4 = ax.bar(x + width, f1_scores["signature_cv"], width,
                    yerr=stds["signature_cv"], capsize=5,
                    color="black", label="Signature (CV)")

    # Test bars
    rects5 = ax.bar(x + width*2, f1_scores["all"], width,
                    color="white", edgecolor="black", label="All (test)")
    rects6 = ax.bar(x + width*3, f1_scores["signature"], width,
                    color="black", label="Signature (test)")
    ax.scatter(x + width*2, dummies["all"], color="red", zorder=3)
    ax.scatter(x + width*3, dummies["signature"], color="red", zorder=3)

    ax.set_xticks(x + width/2)
    ax.set_xticklabels(x_ticks, rotation=45, ha="right")
    ax.set_ylabel("F1 Score")
    ax.legend(ncol=3)

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def plot_heatmap(ordered_taxa,
                 importances,
                 log_fcs,
                 p_values,
                 save_as="../data/plots/heatmap"):
    importances = np.array(importances)
    log_fcs = np.array(log_fcs)
    p_values = np.array(p_values)

    neg_log_p = -np.log10(p_values + 1e-300)  # avoid log(0)

    data = np.vstack([importances, log_fcs, neg_log_p])

    fig, ax = plt.subplots(figsize=(len(ordered_taxa) * 0.4, 3))

    cmap_importance = sns.color_palette("RdBu_r", as_cmap=True)
    cmap_logfc = sns.color_palette("PiYG", as_cmap=True)
    cmap_pval = sns.color_palette("Greys", as_cmap=True)

    # Norms for each row
    norm_importance = None
    norm_logfc = None
    norm_pval = LogNorm(vmin=max(neg_log_p.min(), 1e-2), vmax=neg_log_p.max())

    for i, (row, cmap, norm, label) in enumerate(zip(
        data,
        [cmap_importance, cmap_logfc, cmap_pval],
        [norm_importance, norm_logfc, norm_pval],
        ["Importance", "LogFC", "-log10(p)"]
    )):
        row_reshaped = row[np.newaxis, :]
        ax.imshow(row_reshaped, aspect="auto",
                  cmap=cmap, norm=norm,
                  extent=[0, len(ordered_taxa), 2 - i, 3 - i])

        # Print min/max for this row
        print(f"{label}: min={row.min():.3f}, max={row.max():.3f}")

    ax.set_yticks([2.5, 1.5, 0.5])
    ax.set_yticklabels(["Importance", "LogFC", "-log10(p)"])
    ax.set_xticks(np.arange(len(ordered_taxa)) + 0.5)
    ax.set_xticklabels(ordered_taxa, rotation=90)

    ax.set_xlim(0, len(ordered_taxa))
    ax.set_ylim(0, 3)

    plt.tight_layout()
    plt.savefig(f"{save_as}.svg", format="svg")
    plt.savefig(f"{save_as}.png", format="png", dpi=600)
    plt.close()


def remove_empty_samples(ds, signature_common_features):
    print("Before removing empty samples:", len(ds.taxonomy_counts_df))
    
    mask = (ds.taxonomy_counts_df[signature_common_features].sum(axis=1) > 0)
    ds.taxonomy_counts_df = ds.taxonomy_counts_df.loc[mask]
    ds.metadata_df = ds.metadata_df.loc[mask]

    print("After:", len(ds.taxonomy_counts_df))

    return ds


def compute_feature_importance(rf, ds, features):
    print("Computing feature importance")
    
    result = permutation_importance(rf,
                                    ds.taxonomy_counts_df[features],
                                    ds.metadata_df.loc[ds.taxonomy_counts_df.index, "Treatment"],
                                    n_repeats=25,
                                    random_state=42,
                                    n_jobs=1)
    top_features = {ds.taxonomy_counts_df[features].columns[i]: result.importances_mean[i]
                    for i in range(len(result.importances_mean))}   

    return top_features

def main():
    f1_scores = {"all": [],
                 "signature": [],
                 "all_train": [],
                 "signature_train": [],
                 "all_cv": [],
                 "signature_cv": []}
    dummies = {"all": [],
               "signature": [],
               "all_train": [],
               "signature_train": []}
    stds = {"all_cv": [],
            "signature_cv": []}
    for compartment, study in [("Endosphere", "swift2024drought"),
                               ("Rhizosphere", "zhang2022cross"),
                               ("Rhizosphere", "munoz-ucros2021drought"),
                               ("Bulk Soil", "moore2023microbial")]:
        print("Study:", study)
        
        # Load batch corrected drought dataset for the compartment
        pkl_file = f"../data/merged_l6_{compartment.lower().replace(' ', '_')}_intersection_batch_corrected.pkl"
        with open(pkl_file, 'rb') as handle:
            drought_ds = pkl.load(handle)

        # Load study dataset
        pkl_file = f"../data/{study}_l6_filtered.pkl"
        with open(pkl_file, 'rb') as handle:
            study_ds = pkl.load(handle)

        # Load signature for this compartment
        signature_file = f"../data/maaslin_signature_{compartment.lower().replace(' ', '_')}.pkl"
        with open(signature_file, 'rb') as handle:
            signature = pkl.load(handle)

        # Compute filtered signature (pval <= 0.05)
        filtered_signature = {taxon: signature[taxon]
                              for taxon in signature
                              if signature[taxon]["adj.P.Val"] <= 0.05}

        # Determine feature spaces
        # One feature space: intersection of feature between drought_ds and study_ds
        # The other feature space: intersection between first feature space and the filtered_signature
        drought_features =  set(drought_ds.get_counts_features_columns())
        study_features = set(study_ds.get_counts_features_columns())
        common = list(drought_features.intersection(study_features))
        filtered_signature_common_features = list(set(filtered_signature.keys()).intersection(set(common)))
        
        rf, f1_score, f1_score_train, dummy_score, dummy_score_train, score_cv, cv_std = train_and_evaluate_rf(drought_ds,
                                                                                                               study_ds,
                                                                                                               common)
        f1_scores["all"].append(f1_score)
        f1_scores["all_train"].append(f1_score_train)
        f1_scores["all_cv"].append(score_cv)
        dummies["all"].append(dummy_score)
        dummies["all_train"].append(dummy_score_train)
        stds["all_cv"].append(cv_std)

        # Determine important features based on permutation importance
        top_features = compute_feature_importance(rf, drought_ds, common)
        print(top_features)

        # Train another random forest only with signature taxa
        rf_signature, f1_score, f1_score_train, dummy_score, dummy_score_train, score_cv, cv_std = train_and_evaluate_rf(drought_ds,
                                                                                                                         study_ds,
                                                                                                                         filtered_signature_common_features)
        f1_scores["signature"].append(f1_score)
        f1_scores["signature_train"].append(f1_score_train)
        f1_scores["signature_cv"].append(score_cv)
        dummies["signature"].append(dummy_score)
        dummies["signature_train"].append(dummy_score_train)
        stds["signature_cv"].append(cv_std)
        

        # Intersect with the common features and then
        # get top 25 taxa in the signature, based on p-values
        signature_common_features = list(set(signature.keys()).intersection(set(common)))
        signature = {taxon: signature[taxon]
                     for taxon in signature_common_features}
        sorted_taxa = sorted(signature.items(),
                             key=lambda x: x[1]["adj.P.Val"])
        top25 = dict(sorted_taxa[:25])

        plot_heatmap(top25,
                     [top_features[feature] for feature in top25],
                     [signature[feature]["logFC"] for feature in top25],
                     [signature[feature]["adj.P.Val"] for feature in top25],
                     save_as=f"../data/plots/heatmap_{study}")

    x_ticks = ["swift2024drought",
               "zhang2022cross",
               "munoz-ucros2021drought",
               "moore2023microbial"]
    
    make_grouped_bar_plot(f1_scores,
                          dummies,
                          stds,
                          x_ticks)
                          

if __name__ == "__main__":
    main()
