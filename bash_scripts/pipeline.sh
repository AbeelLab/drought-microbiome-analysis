#!/bin/sh

# Global config variables
data_path="../data"
classifier_path="../silva-138-99-nb-classifier.qza"

# Run for one study
study_id=$1

# Define which steps in the pipeline to run
# Downloading sequencing data and metadata
run_step1=${2:-true}

# Cutadapt to remove primers
run_step2=${3:-false}

# Extract trimming parameters with FIGARO
run_step3=${5:-false}

# DADA2 denoising and ASV assignment
run_step4=${5:-false}

# Taxonomic classification
run_step5=${6:-false}

# Dependency chain if running multiple steps
dependencies=""

update_dependencies() {
    latest_job="$(squeue --me --sort=+i --format="%i" | tail -n 1)"
    
    # remove job array IDs from array jobs
    if [[ "$latest_job" == *"_"* ]]; then
        job_id="$(echo "$latest_job" | cut -d'_' -f1)"
    else
        job_id="${latest_job}"
    fi

    if [ -z "$dependencies" ]; then
        dependencies="afterany:${job_id}"
    else
        dependencies="${dependencies},${job_id}"
    fi

    echo "Dependencies: ${dependencies}"
}

# Run an initial filtering step based on sample metadata
sbatch initial_filter.sbatch $data_path $study_id
update_dependencies

last_batch=1
#$(find "${data_path}/${study_id}/data_batches/" -type f -name "sras_batch*.tsv" -printf "%f\n" | grep -o '[0-9]\+' | sort -n | tail -1)

# Download sequencing data + metadata in batches
if [ $run_step1 = true ] ; then
    sbatch --array="1-${last_batch}" \
	   --dependency=${dependencies} \
	   download.sbatch $data_path \
	                   $study_id

    update_dependencies
fi

# Remove primers if they are present in the reads
if [ $run_step2 = true ] ; then
    sbatch  --array="1-${last_batch}" \
            --dependency=${dependencies} \
	    cutadapt.sbatch $data_path \
                            $study_id 

    update_dependencies
fi

# Run FIGARO to find trimming parameters
if [ $run_step3 = true ] ; then
    # Find trimming parameters per batch
    sbatch  --array="1-${last_batch}" \
            --dependency=${dependencies} \
        figaro.sbatch $data_path \
                      $study_id 

    update_dependencies
fi

# Run DADA2
if [ $run_step4 = true ] ; then
    sbatch  --array="1-${last_batch}" \
            --dependency=${dependencies} \
        dada2.sbatch $data_path \
                     $study_id 

    update_dependencies
fi

# Run taxonomic profiling
if [ $run_step5 = true ] ; then
    sbatch --dependency=${dependencies} \
        classify_taxonomy.sbatch $data_path \
                                 $classifier_path \
                                 $study_id

    update_dependencies
fi
