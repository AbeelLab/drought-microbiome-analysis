#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

# Global arguments
data_path=$1

study_id=$2
trunc_len_f=$3
trunc_len_r=$4
trim_left_f=$5
trim_left_r=$6
max_ee_f=$7
max_ee_r=$8
trunc_q=$9
overlap=${10}

n_threads=${11}

# Output and input directory
working_dir="${data_path}/${study_id}/qiime-dir"

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

qiime dada2 denoise-paired \
      --p-trunc-len-f $trunc_len_f \
      --p-trunc-len-r $trunc_len_r \
      --p-trim-left-f $trim_left_f \
      --p-trim-left-r $trim_left_r \
      --p-max-ee-f $max_ee_f \
      --p-max-ee-r $max_ee_r \
      --p-trunc-q $trunc_q \
      --p-n-threads $n_threads \
      --p-min-overlap $overlap \
      --verbose \
      --i-demultiplexed-seqs "${working_dir}/paired_reads.qza" \
      --o-representative-sequences "${working_dir}/dada2_rep_set.qza" \
      --o-denoising-stats "${working_dir}/dada2_stats.qza" \
      --o-table "${working_dir}/dada2_table.qza"
