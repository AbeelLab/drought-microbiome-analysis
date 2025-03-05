#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --time=00:20:00
#SBATCH --mem=8GB
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

# Global arguments
data_path=$1
classifier_path=$2

# Study-specific input arguments
study_id=$3

# Output and input directory
working_dir="${data_path}/${study_id}/qiime-dir"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Assign taxonomy to ASVs
qiime feature-classifier classify-sklearn \
      --i-classifier $classifier_path \
      --i-reads "${working_dir}/dada2_rep_set.qza" \
      --o-classification "${working_dir}/taxonomy.qza"

# Visualize taxonomy per sample
qiime taxa barplot --i-table "${working_dir}/dada2_table.qza" \
      --i-taxonomy "${working_dir}/taxonomy.qza" \
      --o-visualization "${working_dir}/taxa-bar-plots.qzv"

# Collapse taxonomy per level
# Levels: phylum (2), class (3), order (4), family (5), genus (6)
# Outputs FeatureTable[Frequency], required for ancombc input
for level in $(seq 2 6);
do
    qiime taxa collapse --i-table "${working_dir}/dada2_table.qza" \
          --i-taxonomy "${working_dir}/taxonomy.qza" \
          --p-level $level \
          --o-collapsed-table "${working_dir}/composition_table_l${level}.qza"
    # Extract taxonomy data
    qiime tools extract \
          --input-path "${working_dir}/composition_table_l${level}.qza" \
          --output-path "${working_dir}/composition_table_l${level}"
done
