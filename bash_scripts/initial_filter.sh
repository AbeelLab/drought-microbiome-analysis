#!/bin/sh


data_path=$1
study_id=$2
batch_size=$3

# Output and input directory
working_dir="${data_path}/${study_id}"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Set up QIIME2 cache
cache_dir="/tmp/cache"
qiime tools cache-create --cache $cache_dir
export TMPDIR=$cache_dir

# Download all metadata, extract and filter
# Import accession(s) as qza
qiime tools import \
      --type NCBIAccessionIDs \
      --input-path "${working_dir}/accession.tsv" \
      --output-path "${working_dir}/accession.qza" 

# Download data from NCBI
qiime fondue get-metadata \
      --i-accession-ids "${working_dir}/accession.qza" \
      --p-email bmcosma@tudelft.nl \
      --verbose \
      --o-metadata "${working_dir}/metadata.qza" \
      --o-failed-runs "${working_dir}/failed_IDs.qza" \
      --use-cache $cache_dir

# Visualize metadata
qiime metadata tabulate \
      --m-input-file "${working_dir}/metadata.qza" \
      --o-visualization "${working_dir}/metadata.qzv" \
      --use-cache $cache_dir

# Extract metadata
qiime tools extract \
      --input-path "${working_dir}/metadata.qza" \
      --output-path "${working_dir}/metadata"

conda deactivate
conda activate general

# Filter metadata and write sample names to file
cd ../python_scripts/
python3 initial_filter.py --study $study_id

# Shuffle accessions and split into batches
python3 shuffle_and_batch_sras.py --study $study_id --batch_size $batch_size

cd ../bash_scripts/

rm -rf $cache_dir
