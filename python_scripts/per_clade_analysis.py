import itertools
import os
import pandas as pd
import pickle as pkl

from analysis import run_maaslin_diff_abundance
from copy import deepcopy

def compartment_da_clade_analysis(ds, compartment, clade_pair):
    # Keep only the compartment
    ds.filter_samples_based_on_condition(lambda df: df["RootCompartment"] == compartment)

    # Keep only the two clades
    ds.filter_samples_based_on_condition(lambda df: df["Host"].isin(clade_pair))

    ds.remove_zero_features()

    # If one or both clades don't appear, return
    # This shouldn't happen
    if ds.get_metadata()["Host"].nunique() != 2:
        print(f"Can't compare {compartment} samples for {clade_pair}")
        return

    meta = ds.get_metadata()

    host_counts = meta["Host"].value_counts()
    study_counts = meta.groupby("Host")["StudyID"].nunique()

    print(f"{compartment} | "
          f"{clade_pair[0]}: {host_counts.get(clade_pair[0], 0)} samples "
          f"({study_counts.get(clade_pair[0], 0)} studies), "
          f"{clade_pair[1]}: {host_counts.get(clade_pair[1], 0)} samples "
          f"({study_counts.get(clade_pair[1], 0)} studies)")
    

    # Run DA on batch-corrected dataset
    try:
        results = run_maaslin_diff_abundance(ds,
                                             fixed_effects="Host",
                                             base=clade_pair[0],
                                             condition=clade_pair[1])
        df_results = pd.DataFrame.from_dict(results, orient="index")[["logFC", "adj.P.Val"]]
        df_results = df_results.sort_values("adj.P.Val", ascending=True)

        comparison_label = f"{clade_pair[0]} vs {clade_pair[1]} (reference: {clade_pair[0]})"
        filename = f"DA_{compartment.lower()}_{clade_pair[0]}_vs_{clade_pair[1]}_ref_{clade_pair[0].replace(' ', '_')}.csv"
        filepath = os.path.join("../raw_data", filename)

        with open(filepath, "w") as f:
            with open(filepath, "w") as f:
                f.write(f"# {comparison_label}\n")
            df_results.to_csv(filepath, mode="a", header=True)

        print(f"Saved: {filepath}")
    except:
        print(f"Failed: {clade_pair}")
    

def main():
    clades = ["BOP clade", "PACMAD clade", "Asterids", "Fabids", "Malvids"]
    clade_pairs = list(itertools.combinations(clades, 2))

    # DA analysis between clades, per compartment
    pkl_file = f"../datasets/phylum_batch_corrected.pkl"
    with open(pkl_file, 'rb') as handle:
        batch_corrected_phyla = pkl.load(handle)
    for compartment in ["Endosphere", "Rhizosphere"]:
        print("---", compartment)
        for clade_pair in clade_pairs:
            compartment_da_clade_analysis(ds=deepcopy(batch_corrected_phyla),
                                          compartment=compartment,
                                          clade_pair=clade_pair)


    clades = ["BOP clade", "PACMAD clade", "Asterids", "Fabids", "Malvids"]
    thresholds = [0.75, 0.8, 0.85, 0.9, 0.95]

    for compartment in ["Endosphere", "Rhizosphere"]:
        pkl_file = f"../datasets/{compartment.lower().replace(' ', '_')}/genus.pkl"
        with open(pkl_file, 'rb') as handle:
            ds = pkl.load(handle)

        # Store results per clade
        results = []

        for clade in clades:
            ds_copy = deepcopy(ds)
            ds_copy.filter_samples_based_on_condition(lambda df: df["Host"] == clade)

            num_samples = len(ds_copy.get_counts_features())
            core_sizes = []

            for threshold in thresholds:
                core, _ = ds_copy.get_core(min_prevalence=threshold)
                core_sizes.append(len(core))

            results.append([num_samples] + core_sizes)

        columns = ["num_samples"] + [f"core_{int(th*100)}" for th in thresholds]
        df_results = pd.DataFrame(results, index=clades, columns=columns)

        filename = f"clade_cores_{compartment.lower()}.csv"
        filepath = os.path.join("../raw_data", filename)
        df_results.to_csv(filepath)

        print(f"Saved core sizes for {compartment} to {filepath}")
            

if __name__ == "__main__":
    main()
