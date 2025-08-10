#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out
#SBATCH --qos=medium
#SBATCH --time=24:00:00
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
# Final trimming parameters should be logged
trim_param_file="${study_dir}/figaro_trim_params.txt"
# For these two studies FIGARO does terminate, but the resulting taxonomic profiling is not good
if [[ -f $trim_param_file && \
      "$study_id" != "simmons2020drought" && \
      "$study_id" != "azarbad2020four" && \
      "$study_id" != "swift2024drought" && \
      "$study_id" != "faist2022potato" && \
      "$study_id" != "azarbad2022response" ]]; then
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


# Defaults
forward_trim=0
reverse_trim=0
trunc_q=2
forward_max_ee=2
reverse_max_ee=2
overlap=12

# For this particular study neither figaro parameters nor the 90% amplicon length work
if [[ "$study_id" == "moore2023microbial" ]]; then
    forward_trunc=240
    reverse_trunc=200
    forward_trim=25
    reverse_trim=26
    trunc_q=2
    forward_max_ee=2
    reverse_max_ee=4
    overlap=100
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
     --p-trim-left-f $forward_trim \
     --p-trim-left-r $reverse_trim \
     --p-max-ee-f $forward_max_ee \
     --p-max-ee-r $reverse_max_ee \
     --p-trunc-q $trunc_q \
     --p-min-overlap $overlap \
     --p-n-threads 32 \
     --verbose \
     --p-no-hashed-feature-ids \
     --i-demultiplexed-seqs "${qiime_dir}/paired_reads_no_primers.qza" \
     --o-representative-sequences "${qiime_dir}/dada2_rep_set.qza" \
     --o-denoising-stats "${qiime_dir}/dada2_stats.qza" \
     --o-table "${qiime_dir}/dada2_table.qza" \
     --use-cache $cache_dir
