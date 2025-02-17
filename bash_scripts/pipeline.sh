#!/bin/sh

# Global config variables
global_path="/tudelft.net/staff-umbrella/abeellab/bmcosma/resilience-characterization"
data_path="${global_path}/data"
classifier_path="${global_path}/silva-138-99-nb-classifier.qza"

# DADA2 defaults
trim_left_f=0
trim_left_r=0
max_ee_f="2.0"
max_ee_r="2.0"
trunc_q=2
overlap=12

# Run for one study
study_id=$1

# Load study config
# DADA2 defaults should be overwritten if specified in study config
source "${study_id}_config"

# Define which steps in the pipeline to run
# By default: yes
run_step1=${2:-true}
run_step2=${3:-true}

# Dependency chain if running multiple steps
dependencies=""

update_dependencies() {
    sleep 15s
    #https://stackoverflow.com/questions/73773769/slurm-get-job-id-of-last-run-jobs
    latest_job="$(squeue --me --sort=+i --format="%i" | tail -n 1)"
    if [ -z "$dependencies" ]; then
        dependencies="afterok:${latest_job}"
    else
        dependencies="${dependencies},${latest_job}"
    fi
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

echo "Dependencies: ${dependencies}"

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

echo "Dependencies: ${dependencies}"

