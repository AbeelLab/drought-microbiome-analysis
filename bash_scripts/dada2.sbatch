#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out
#SBATCH --qos=short
#SBATCH --time=04:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=32
#SBATCH --gres=jobspernode:one:1
#SBATCH --tmp=150G

data_path=$1
study_id=$2

# Main output and input directory
study_dir="${data_path}/${study_id}"
working_dir="${study_dir}/data_batches"
# Batch-specific directories and files
qiime_dir="${working_dir}/qiime-dir-batch${SLURM_ARRAY_TASK_ID}"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh

# Determine trimming parameters as medians of batches
conda activate general

cd ../python_scripts/
python3 extract_figaro_parameters.py --study $study_id

cd ../bash_scripts/

conda deactivate

# Load trimming parameters
trim_param_file="${study_dir}/figaro_trim_params.txt"
forward_trunc=$(sed -n '1p' "$trim_param_file")
reverse_trunc=$(sed -n '2p' "$trim_param_file")

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

conda activate qiime2-amplicon-2024.10

# Set up QIIME2 cache
cache_dir="/tmp/job_${SLURM_JOB_ID}"
srun qiime tools cache-create --cache $cache_dir
export TMPDIR=$cache_dir

# Clean up cache after job execution
function clean_up {
    echo "Before cleanup:"
    find /tmp -user bmcosma -print0 | du -chs --files0-from=-
    echo "---"
    
    rm -rf $cache_dir

    echo "After cleanup:"
    find /tmp -user bmcosma -print0 | du -chs --files0-from=-
    exit
}

trap 'clean_up' EXIT

# Do not hash feature IDs to ensure the same ASVs from different batches
# will be merged correctly.
# https://forum.qiime2.org/t/merging-dada2-feature-tables-what-about-likely-duplicate-features-in-different-tables/12419/2
srun qiime dada2 denoise-paired \
     --p-trunc-len-f $forward_trunc \
     --p-trunc-len-r $reverse_trunc \
     --p-n-threads 32 \
     --verbose \
     --p-no-hashed-feature-ids \
     --i-demultiplexed-seqs "${qiime_dir}/paired_reads_no_primers.qza" \
     --o-representative-sequences "${qiime_dir}/dada2_rep_set.qza" \
     --o-denoising-stats "${qiime_dir}/dada2_stats.qza" \
     --o-table "${qiime_dir}/dada2_table.qza" \
     --use-cache $cache_dir
