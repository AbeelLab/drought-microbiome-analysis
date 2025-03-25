#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out
#SBATCH --qos=short
#SBATCH --time=00:15:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=8

# Global arguments
data_path=$1
study_id=$2

# Main output and input directory
working_dir="${data_path}/${study_id}/data_batches"
# Batch-specific directories and files
qiime_dir="${working_dir}/qiime-dir-batch${SLURM_ARRAY_TASK_ID}"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Extract QIIME2 reads artifact
qiime tools extract \
      --input-path "${qiime_dir}/paired_reads_no_primers.qza" \
      --output-path "${qiime_dir}/paired_reads_no_primers"

conda deactivate
conda activate general

config_file="../python_scripts/config.yml"
amplicon_length=$(python -c "import yaml; print(yaml.safe_load(open('$config_file'))['$study_id']['amplicon_length'])")

conda deactivate
conda activate figaro

input_dir=$(echo "${qiime_dir}/paired_reads_no_primers/"*/data)

function clean_up {
    rm -rf $input_dir
    exit
}

trap 'clean_up' EXIT

srun python3 ../figaro/figaro/figaro.py -i $input_dir \
     -o "${qiime_dir}/figaro-output" \
     -a $amplicon_length \
     -f 0 \
     -r 0 \
     -F illumina
