#!/bin/sh

# Global config variables
data_path="../data"
classifier_path="../silva-138-99-nb-classifier.qza"

# DADA2 conservative parameters
trunc_len_f=220
trunc_len_r=170
trim_left_f=25
trim_left_r=26
max_ee_f="2.0"
max_ee_r="4.0"
trunc_q=2
overlap=100

# Run for one study
study_id=$1

# Load study config
# Will overwirte DADA2 parameters where specified in the config
source "configs/${study_id}_config"

# Define which steps in the pipeline to run
# By default: all true

# Downloading sequencing data and metadata
run_step1=${2:-true}

# DADA2 denoising and ASV assignment
run_step2=${3:-true}

# Taxonomic classification
run_step3=${4:-true}

# Dependency chain if running multiple steps
dependencies=""

update_dependencies() {
    sleep 3s
    latest_job="$(squeue --me --sort=+i --format="%i" | tail -n 1)"
    if [ -z "$dependencies" ]; then
        dependencies="afterok:${latest_job}"
    else
        dependencies="${dependencies},${latest_job}"
    fi

    echo "Dependencies: ${dependencies}"
}

# Download data
if [ $run_step1 = true ] ; then
    sbatch  --qos=$qos_download \
            --time=$time_download \
    download.sbatch $data_path \
                    $study_id \
                    $accession

    update_dependencies
fi

# Run DADA2
if [ $run_step2 = true ] ; then
    sbatch  --qos=$qos_dada2 \
            --time=$time_dada2 \
            --cpus-per-task=$n_threads_dada2  \
            --mem=$mem_dada2 \
            --dependency=${dependencies} \
        dada2.sbatch $data_path \
                     $study_id \
                     $trunc_len_f \
                     $trunc_len_r \
                     $trim_left_f \
                     $trim_left_r \
                     $max_ee_f \
                     $max_ee_r \
                     $trunc_q \
                     $overlap \
                     $n_threads_dada2

    update_dependencies
fi

# Run taxonomic profiling
if [ $run_step3 = true ] ; then
    sbatch --dependency=${dependencies} \
        classify_taxonomy.sbatch $data_path \
                                 $classifier_path \
                                 $study_id

    update_dependencies
fi
