#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out
#SBATCH --qos=medium
#SBATCH --time=08:00:00

data_path=$1
study_id=$2

# Main output and input directory
working_dir="${data_path}/${study_id}/data_batches"
# Batch-specific directories and files
tsv_sra_list="${working_dir}/sras_batch${SLURM_ARRAY_TASK_ID}.tsv"
qza_sra_list="${working_dir}/sras_batch${SLURM_ARRAY_TASK_ID}.qza"
qiime_dir="${working_dir}/qiime-dir-batch${SLURM_ARRAY_TASK_ID}"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Set up QIIME2 cache
cache_dir="../tmp/job_${SLURM_JOB_ID}"
qiime tools cache-create --cache $cache_dir
# Clean up cache after job execution
trap 'rm -rf "$cache_dir"' EXIT

# Import accession(s) as qza
qiime tools import \
      --type NCBIAccessionIDs \
      --input-path $tsv_sra_list \
      --output-path $qza_sra_list

# Remove qiime-dir if it already exists (otherwise q2-fondue doesn't work)
rm -rf $qiime_dir

# Download data from NCBI
qiime fondue get-all \
      --i-accession-ids $qza_sra_list \
      --p-email bmcosma@tudelft.nl \
      --p-retries 5 \
      --verbose \
      --use-cache $cache_dir \
      --output-dir $qiime_dir

# Extract metadata for later processing
qiime tools extract \
      --input-path "${qiime_dir}/metadata.qza" \
      --output-path "${qiime_dir}/metadata"

# Visualize metadata
qiime metadata tabulate \
      --m-input-file "${qiime_dir}/metadata.qza" \
      --use-cache $cache_dir \
      --o-visualization "${qiime_dir}/metadata.qzv" 

# Visualize read statistics
qiime demux summarize \
      --i-data "${qiime_dir}/paired_reads.qza" \
      --use-cache $cache_dir \
      --o-visualization "${qiime_dir}/paired_reads.qzv"
