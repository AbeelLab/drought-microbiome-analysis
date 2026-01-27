import pandas as pd
import pickle as pkl

from copy import deepcopy
from skbio.diversity import alpha_diversity
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

def compute_shannon(ds, treatment, clade):
    if treatment != "All":
        ds.filter_samples_based_on_condition(lambda df: df["Treatment"] == treatment)
    if clade != "All" and clade is not None:
        ds.filter_samples_based_on_condition(lambda df: df["Host"] == clade)

    df = ds.get_counts_features()
    shannon = alpha_diversity(metric="shannon",
                              counts=df.values,
                              ids=df.index)

    return shannon
        
    
def main():
    compartments = ["Endosphere", "Rhizosphere", "Bulk soil"]
    clades = ["All", "BOP clade", "PACMAD clade", "Asterids", "Fabids", "Malvids"]
    treatments = ["All", "Drought", "Control"]

    results = []
    pvals = []

    for compartment in compartments:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        # We should batch correct.
        # MMUPHin will not remove zeros (required for Shannon diversity)
        # but batch correction should be done for rank-sum test.
        ds.apply_mmuphin_be_correction()

        if compartment != "Bulk soil":
            clade_iter = clades
        else:
            clade_iter = [None]

        for clade in clade_iter:
            label = (compartment if clade in ["All", None]
                     else f"{compartment} | {clade}")

            sh_all = compute_shannon(deepcopy(ds), "All", clade)
            sh_drought = compute_shannon(deepcopy(ds), "Drought", clade)
            sh_control = compute_shannon(deepcopy(ds), "Control", clade)

            row = {"Sample type": label,
                   "All": sh_all.mean(),
                   "Only drought": sh_drought.mean(),
                   "Only control": sh_control.mean(),
                   "n_drought": len(sh_drought),
                   "n_control": len(sh_control)}

            stat, pval = mannwhitneyu(sh_drought.values,
                                      sh_control.values,
                                      alternative="two-sided")

            row["p_value"] = pval
            pvals.append(pval)

            results.append(row)

    df_results = pd.DataFrame(results)

    df_results["p_adj"] = multipletests(df_results["p_value"],
                                        method="fdr_bh")[1]
    
    df_results.to_csv("../raw_data/shannon_results.csv", index=False)


if __name__ == "__main__":
    main()
