#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --time=00:15:00
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

data_path=$1
classifier_path=$2
study_id=$3

# Main output and input directory
study_dir="${data_path}/${study_id}"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Merge dada2 files across batches
dada2_rep_files=()
for batch in "${study_dir}/data_batches/qiime-dir-batch"*; do
    if [ -f "${batch}/dada2_rep_set.qza" ]; then
        dada2_rep_files+=( "${batch}/dada2_rep_set.qza" )
    fi
done

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

srun qiime feature-table merge \
     --i-tables "${dada2_rep_files[@]}" \
     --p-overlap-method sum \
     --o-merged-table "${study_dir}/dada2_rep_set.qza" \
     --use-cache $cache_dir

# Assign taxonomy to ASVs
srun qiime feature-classifier classify-sklearn \
     --i-classifier $classifier_path \
     --i-reads "${study_dir}/dada2_rep_set.qza" \
     --o-classification "${study_dir}/taxonomy.qza" \
     --use-cache $cache_dir

# Visualize taxonomy per sample
srun qiime taxa barplot --i-table "${study_dir}/dada2_table.qza" \
     --i-taxonomy "${study_dir}/taxonomy.qza" \
     --o-visualization "${study_dir}/taxa-bar-plots.qzv" \
     --use-cache $cache_dir

# Collapse taxonomy per level
# Levels: phylum (2), class (3), order (4), family (5), genus (6)
# Outputs FeatureTable[Frequency], required for ancombc input
for level in $(seq 2 6);
do
    srun qiime taxa collapse --i-table "${study_dir}/dada2_table.qza" \
         --i-taxonomy "${study_dir}/taxonomy.qza" \
         --p-level $level \
         --o-collapsed-table "${study_dir}/composition_table_l${level}.qza" \
	 --use-cache $cache_dir
    # Extract taxonomy data
    srun qiime tools extract \
         --input-path "${study_dir}/composition_table_l${level}.qza" \
         --output-path "${study_dir}/composition_table_l${level}"
done
