#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

# Global arguments
data_path=$1

# Study-specific input arguments
study_id=$2
accession=$3

# Output and input directory
working_dir="${data_path}/${study_id}"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Import accession(s) as qza
qiime tools import \
      --type NCBIAccessionIDs \
      --input-path "${working_dir}/${accession}.tsv" \
      --output-path "${working_dir}/${accession}.qza"

# Download data from NCBI
qiime fondue get-all \
      --i-accession-ids "${working_dir}/${accession}.qza" \
      --p-email bmcosma@tudelft.nl \
      --p-retries 5 \
      --verbose \
      --output-dir "${working_dir}/qiime-dir"

# Visualize metadata
qiime metadata tabulate \
      --m-input-file "${working_dir}/qiime-dir/metadata.qza" \
      --o-visualization "${working_dir}/qiime-dir/metadata.qzv"

# Visualize read statistics
qiime demux summarize \
      --i-data "${working_dir}/qiime-dir/paired_reads.qza" \
      --o-visualization "${working_dir}/qiime-dir/paired_reads.qzv"
