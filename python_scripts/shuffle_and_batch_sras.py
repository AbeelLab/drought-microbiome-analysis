import yaml
import os
import argparse
import random
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study")
    parser.add_argument("--batch_size")
    args = parser.parse_args()

    study = args.study
    batch_size = int(args.batch_size)

    # Set random seed for shuffling
    random.seed(42)
    
    with open('config.yml') as f:
        config = yaml.safe_load(f)

    study_path = os.path.join(config["data_path"], study)
    accession_file = os.path.join(study_path, "filtered_sras.tsv")

    with open(accession_file, "r") as f:
        lines = f.read().splitlines()

    if not lines:
        print(f"No data found in {accession_file}")
        return
        
    # accession.tsv:
    # id
    # SRA...
    # SRA...
    # ...
    # SRA...
    header = lines[0]
    sra_ids = lines[1:]

    # Shuffle SRA IDs
    # If we process data in batches they should be randomly distributed
    random.shuffle(sra_ids)

    shuffled_file = os.path.join(study_path, "sras_shuffled.tsv")
    with open(shuffled_file, "w") as f:
        f.write(header + "\n")
        for sra in sra_ids:
            f.write(sra + "\n")
    print(f"Shuffled accession file saved as: {shuffled_file}")

    batches_dir = os.path.join(study_path, "data_batches")
    # remove batches from previous runs
    if os.path.exists(batches_dir):
        shutil.rmtree(batches_dir)
    os.makedirs(batches_dir, exist_ok=True)

    batches = []
    for i in range(0, len(sra_ids), batch_size):
        batch = sra_ids[i:i+batch_size]
        batches.append(batch)

    # If last batch is less than 15 samples, merge with second-to-last
    if len(batches) > 1 and len(batches[-1]) < 15:
        batches[-2].extend(batches[-1])
        batches.pop()
    
    for idx, batch in enumerate(batches, start=1):
        batch_file = os.path.join(batches_dir, f"sras_batch{idx}.tsv")
        with open(batch_file, "w") as f:
            f.write(header + "\n")
            for sra in batch:
                f.write(sra + "\n")
        print(f"Batch file created: {batch_file}")

if __name__ == "__main__":
    main()
