import pandas as pd
import pickle
import sys

# Specify the desired row order
ordered_keys = [
    "xu2018drought",
    "simmons2020drought",
    "fitzpatrick2018assembly",
    "azarbad2020four",
    "santos-medellin2021prolonged",
    "hoefle2024oak",
    "bandopadhyay2024disentangling",
    "naylor2017drought",
    "faist2022potato",
    "santos-medellin2017drought",
    "azarbad2022response",
    "wipf2021distinguishing",
    "moore2023microbial",
    "swift2024drought",
    "zhang2022cross",
    "munoz-ucros2021drought"
]

def main(pickle_path):
    # Load the dictionary from pickle
    with open(pickle_path, 'rb') as f:
        data = pickle.load(f)

    # Create DataFrame
    df = pd.DataFrame.from_dict(data, orient='index')

    # Reorder the rows according to ordered_keys
    # Filter to keep only those in the dict (in case some are missing)
    present_keys = [key for key in ordered_keys if key in df.index]
    df_ordered = df.loc[present_keys]

    # Save to CSV
    output_csv = pickle_path.replace('.pkl', '.csv')
    df_ordered.to_csv(output_csv)
    print(f"Saved CSV to: {output_csv}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python convert_log_to_csv.py <samples_log.pkl | features_log.pkl>")
        sys.exit(1)

    main(sys.argv[1])
