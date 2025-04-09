#!/bin/sh

data_path=$1
study_id=$2
# One accession or comma-separated: accession_1,accession_2,...
accession=$3

dir="${data_path}/${study_id}"
mkdir -p "$dir"

file="${dir}/accession.tsv"

# Header (required by q2-fondue)
echo "id" > "$file"

# Split accessions on commas and write on separate rows
echo "$accession" | tr ',' '\n' >> "$file"
