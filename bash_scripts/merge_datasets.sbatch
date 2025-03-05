#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --time=00:05:00
#SBATCH --mem=4GB
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

# Global arguments
data_path=$1

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate general

# Merge per level
for level in $(seq 2 6);
do    
    python3 ../python_scripts/merge_datasets.py $data_path $level
done

