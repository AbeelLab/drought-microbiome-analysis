import json
import os
import yaml
import argparse
import glob
import statistics

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study")
    args = parser.parse_args()

    study = args.study

    with open('config.yml') as f:
        config = yaml.safe_load(f)
    
    study_dir = os.path.join(config["data_path"],
                             study)
    trim_positions_forward = []
    trim_positions_reverse = []
    
    pattern = os.path.join(study_dir, "data_batches", "qiime-dir-batch*")
    for batch_dir in glob.glob(pattern):
        json_file = os.path.join(batch_dir, "trimParameters.json")

        with open(json_file) as jf:
            data = json.load(jf)
            tp = data[0]["trimPosition"]
            trim_positions_forward.append(tp[0])
            trim_positions_reverse.append(tp[1])
            
    median_forward = statistics.median(trim_positions_forward)
    median_reverse = statistics.median(trim_positions_reverse)

    print("Standard deviation forward:", statistics.stdev(trim_positions_forward))
    print("Standard deviation reverse:", statistics.stdev(trim_positions_reverse))
    
    output_file = os.path.join(study_dir, "figaro_trim_params.txt")
    with open(output_file, "w") as outf:
        outf.write(f"{median_forward}\n{median_reverse}\n")

if __name__ == "__main__":
    main()
