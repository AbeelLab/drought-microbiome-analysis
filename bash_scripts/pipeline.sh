#!/bin/sh

. ./parser.sh "$@"

# If accession is provided, create study directory and .tsv file within this directory required by q2-fondue
if [ -n "$accession" ]; then
    echo "Formatting accession(s) $accession for study $study_id..."
    bash create_study_dir.sh $data_path $study_id $accession
fi

# Keep track of dependency chain
dependencies=""

update_dependencies() {
    latest_job="$(squeue --me --sort=+i --format="%i" | tail -n 1)"
    
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

# TO DO: separate filtering and batching
# And make them tweakable
if [ "$run_initial_filter_and_batching" = true ] ; then
    bash initial_filter.sh "$data_path" "$study_id" "$batch_size"
fi

last_batch=$(find "${data_path}/${study_id}/data_batches/" -type f -name "sras_batch*.tsv" -printf "%f\n" | grep -o '[0-9]\+' | sort -n | tail -1)

if [ "$run_download" = true ] ; then
    echo "Number of batches: ${last_batch}"

    # Create list of batches that have not yet been downloaded
    array_job_list=""
    for i in $(seq 1 $last_batch);
    do
	working_dir="${data_path}/${study_id}/data_batches"
	reads_file="${working_dir}/qiime-dir-batch${i}/paired_reads.qza"

	if [ ! -f $reads_file ]; then
	    array_job_list="${array_job_list},${i}"
	fi
    done

    # Remove first comma
    array_job_list="${array_job_list:1}"

    sbatch --array=$array_job_list \
           download.sbatch \
           "$data_path" \
           "$study_id"
    update_dependencies
fi


if [ "$run_cutadapt" = true ] ; then
    sbatch --array="1-${last_batch}" \
           --dependency=${dependencies} \
           cutadapt.sbatch \
           "$data_path" \
           "$study_id"
    update_dependencies
fi


if [ "$run_figaro" = true ] ; then
    sbatch --array="1-${last_batch}" \
           --dependency=${dependencies} \
           figaro.sbatch \
           "$data_path" \
           "$study_id"
    update_dependencies
fi


if [ "$run_dada2" = true ] ; then
    sbatch --dependency=${dependencies} \
           load_figaro_trim_params.sbatch $data_path $study_id $last_batch

    update_dependencies

    sbatch --array="1-${last_batch}" \
           --dependency=${dependencies} \
           dada2.sbatch \
           "$data_path" \
           "$study_id"
    update_dependencies
fi


if [ "$run_taxonomic_profiling" = true ] ; then
    sbatch --dependency=${dependencies} \
           classify_taxonomy.sbatch \
           "$data_path" \
           "$classifier_path" \
           "$study_id"
fi
