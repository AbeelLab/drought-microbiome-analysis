#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --time=00:30:00
#SBATCH --mem=16GB
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
dada2_tables=()
dada2_rep_files=()
for batch in "${study_dir}/data_batches/qiime-dir-batch"*; do
    dada2_tables+=( "${batch}/dada2_table.qza" )
    dada2_rep_files+=( "${batch}/dada2_rep_set.qza" )
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

# # ---> ASV table merging
srun qiime feature-table merge \
     --i-tables "${dada2_tables[@]}" \
     --p-overlap-method sum \
     --o-merged-table "${study_dir}/dada2_table.qza" \
     --use-cache $cache_dir

# Extract dada2_table
srun qiime tools extract \
     --input-path "${study_dir}/dada2_table.qza" \
     --output-path "${study_dir}/dada2_table"

qiime feature-table merge-seqs \
      --i-data "${dada2_rep_files[@]}" \
      --o-merged-data "${study_dir}/dada2_rep_set.qza" \
      --use-cache $cache_dir


# ---> Taxonomic classification
srun qiime feature-classifier classify-sklearn \
     --i-classifier $classifier_path \
     --i-reads "${study_dir}/dada2_rep_set.qza" \
     --o-classification "${study_dir}/taxonomy.qza" \
     --use-cache $cache_dir


# ---> Fitering
# Filter ASVs based on taxonomy
# Keep only ASVs assigned at least phylum level (p__)
qiime taxa filter-table \
      --i-table "${study_dir}/dada2_table.qza" \
      --i-taxonomy "${study_dir}/taxonomy.qza" \
      --p-mode contains \
      --p-include "p__" \
      --p-exclude "Chloroplast,Mitochondria,Eukaryota,Unassigned,Unclassified" \
      --o-filtered-table "${study_dir}/dada2_asv_filtered_table.qza"

# Filter samples
# Based on number of ASVs
qiime feature-table filter-samples \
      --i-table "${study_dir}/dada2_asv_filtered_table.qza" \
      --p-min-frequency 1000 \
      --o-filtered-table "${study_dir}/dada2_sample_filtered_table.qza"
# Based on number of taxonomic features
qiime feature-table filter-samples \
      --i-table "${study_dir}/dada2_sample_filtered_table.qza" \
      --p-min-features 10 \
      --o-filtered-table "${study_dir}/dada2_final_filtered_table.qza"

# Filter representative sequences
qiime feature-table filter-seqs \
      --i-data "${study_dir}/dada2_rep_set.qza" \
      --i-table "${study_dir}/dada2_final_filtered_table.qza" \
      --o-filtered-data "${study_dir}/dada2_filtered_rep_set.qza"

# Extract tables
srun qiime tools extract \
     --input-path "${study_dir}/dada2_asv_filtered_table.qza" \
     --output-path "${study_dir}/dada2_asv_filtered_table"

srun qiime tools extract \
     --input-path "${study_dir}/dada2_final_filtered_table.qza" \
     --output-path "${study_dir}/dada2_final_filtered_table"

# ---> Visualization
# Visualize taxonomy per sample
srun qiime taxa barplot --i-table "${study_dir}/dada2_final_filtered_table.qza" \
     --i-taxonomy "${study_dir}/taxonomy.qza" \
     --o-visualization "${study_dir}/taxa-bar-plots.qzv" \
     --use-cache $cache_dir

# Collapse taxonomy per level
# Levels: phylum (2), class (3), order (4), family (5), genus (6)
# Outputs FeatureTable[Frequency], required for ancombc input
for level in $(seq 2 6);
do
    srun qiime taxa collapse --i-table "${study_dir}/dada2_final_filtered_table.qza" \
         --i-taxonomy "${study_dir}/taxonomy.qza" \
         --p-level $level \
         --o-collapsed-table "${study_dir}/composition_table_l${level}.qza" \
	 --use-cache $cache_dir

    # To avoid multiple extracted tables, remove existing directories
    rm -rf "${study_dir}/composition_table_l${level}"
    
    # Extract taxonomy data
    srun qiime tools extract \
         --input-path "${study_dir}/composition_table_l${level}.qza" \
         --output-path "${study_dir}/composition_table_l${level}"
done
