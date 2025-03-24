#!/bin/sh
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out
#SBATCH --qos=short
#SBATCH --time=00:30:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=8
#SBATCH --gres=jobspernode:one:1

data_path=$1
study_id=$2

# Main output and input directory
working_dir="${data_path}/${study_id}/data_batches"
# Batch-specific directories and files
qiime_dir="${working_dir}/qiime-dir-batch${SLURM_ARRAY_TASK_ID}"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Set up QIIME2 cache
cache_dir="/tmp/job_${SLURM_JOB_ID}"
export TMPDIR=$cache_dir

qiime tools cache-create --cache $cache_dir

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

qiime cutadapt trim-paired \
      --p-cores 8 \
      --use-cache $cache_dir \
      --i-demultiplexed-sequences "${qiime_dir}/paired_reads.qza" \
      --o-trimmed-sequences "${qiime_dir}/paired_reads_no_primers.qza" \
      --verbose
