#!/bin/sh

data_path="../data"
classifier_path="../silva-138-99-nb-classifier.qza"
batch_size=50

study_id=""
accession=""

run_download=""
run_cutadapt=""
run_dada2=""
run_taxonomic_profiling=""
run_initial_filter_and_batching=""

pipeline_flags_specified=0

while [ "$#" -gt 0 ]; do
    case "$1" in
        --study_id)
            study_id="$2"
            shift 2
            ;;
        --accession)
            accession="$2"
            shift 2
            ;;
        --data_path)
            data_path="$2"
            shift 2
            ;;
        --classifier_path)
            classifier_path="$2"
            shift 2
            ;;
        --batch_size)
            batch_size="$2"
            shift 2
            ;;
        --run_initial_filter_and_batching)
            run_initial_filter_and_batching=true
            pipeline_flags_specified=1
            shift 1
            ;;
        --run_download)
            run_download=true
            pipeline_flags_specified=1
            shift 1
            ;;
        --run_cutadapt)
            run_cutadapt=true
            pipeline_flags_specified=1
            shift 1
            ;;
        --run_figaro)
            run_figaro=true
            pipeline_flags_specified=1
            shift 1
            ;;
        --run_dada2)
            run_dada2=true
            pipeline_flags_specified=1
            shift 1
            ;;
        --run_taxonomic_profiling)
            run_taxonomic_profiling=true
            pipeline_flags_specified=1
            shift 1
            ;;
        --help)
            ./usage.sh
            exit 1
            ;;
        *)
            echo "Unknown parameter: $1"
            ./usage.sh
            exit 1
            ;;
    esac
done

if [ -z "$study_id" ]; then
    echo "Error: --study_id is required"
    ./usage.sh
    exit 1
fi

if [ "$pipeline_flags_specified" -eq 0 ]; then
    run_initial_filter_and_batching=true
    run_download=true
    run_cutadapt=true
    run_dada2=true
    run_taxonomic_profiling=true
fi

export data_path
export classifier_path
export batch_size
export study_id
export accession
export run_initial_filter_and_batching
export run_download
export run_cutadapt
export run_dada2
export run_taxonomic_profiling
