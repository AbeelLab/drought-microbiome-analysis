#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out
#SBATCH --qos=medium
#SBATCH --time=06:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=32
#SBATCH --gres=jobspernode:one:1
#SBATCH --tmp=200G

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

data_path=$1
study_id=$2

# Main output and input directory
study_dir="${data_path}/${study_id}"
working_dir="${study_dir}/data_batches"
# Batch-specific directories and files
qiime_dir="${working_dir}/qiime-dir-batch${SLURM_ARRAY_TASK_ID}"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh

# Load trimming parameters
# TO DO: these if-statements need to be changed and replaced with more consistent rules
# Final trimming parameters should be logged
trim_param_file="${study_dir}/figaro_trim_params.txt"
# For these two studies FIGARO does terminate, but the resulting taxonomic profiling is not good
if [[ -f $trim_param_file &&
	  "$study_id" != "simmons2020drought" &&
	  "$study_id" != "azarbad2020four" ]]; then
    # If FIGARO ran successfully, load optimal parameters
    forward_trunc=$(sed -n '1p' "$trim_param_file")
    reverse_trunc=$(sed -n '2p' "$trim_param_file")
else
    # Otherwise trim 10%
    config_file="../python_scripts/config.yml"
    conda activate general
    amplicon_length=$(python -c "import yaml; print(yaml.safe_load(open('$config_file'))['$study_id']['amplicon_length'])")
    conda deactivate
    forward_trunc=$((amplicon_length - amplicon_length / 10))
    reverse_trunc=$((amplicon_length - amplicon_length / 10))
fi

echo $forward_trunc
echo $reverse_trunc

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
