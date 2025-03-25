#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --time=00:10:00
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

# Global arguments
data_path=$1

# Study-specific input arguments
study_id=$2

# Output and input directory
working_dir="${data_path}/${study_id}"

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Set up QIIME2 cache
cache_dir="/tmp/job_${SLURM_JOB_ID}"
srun qiime tools cache-create --cache $cache_dir
export TMPDIR=$cache_dir

# Clean up cache after job execution
function clean_up {
    echo "Before cleanup:"
    find /tmp -user bmcosma -print0 | du -chs --files0-from=-
    echo "---"
    
    rm -rf $cache_dir

    echo "After cleanup:"
    find /tmp -user bmcosma -print0 | du -chs --files0-from=-
    exit
}

trap 'clean_up' EXIT

# Download all metadata, extract and filter
# Import accession(s) as qza
srun qiime tools import \
     --type NCBIAccessionIDs \
     --input-path "${working_dir}/accession.tsv" \
     --output-path "${working_dir}/accession.qza" 

# Download data from NCBI
srun qiime fondue get-metadata \
     --i-accession-ids "${working_dir}/accession.qza" \
     --p-email bmcosma@tudelft.nl \
     --verbose \
     --o-metadata "${working_dir}/metadata.qza" \
     --o-failed-runs "${working_dir}/failed_IDs.qza" \
     --use-cache $cache_dir

# Visualize metadata
srun qiime metadata tabulate \
     --m-input-file "${working_dir}/metadata.qza" \
     --o-visualization "${working_dir}/metadata.qzv" \
     --use-cache $cache_dir

# Extract metadata
srun qiime tools extract \
     --input-path "${working_dir}/metadata.qza" \
     --output-path "${working_dir}/metadata"

conda deactivate
conda activate general

# Filter metadata and write sample names to file
cd ../python_scripts/
python3 initial_filter.py --study $study_id

# Shuffle accessions and split into batches
python3 shuffle_and_batch_sras.py --study $study_id

cd ../bash_scripts/
