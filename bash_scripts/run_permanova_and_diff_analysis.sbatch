#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --qos=short
#SBATCH --time=04:00:00
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

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

for level in 6; do
    for suffix in "" "_batch_corrected"; do
        echo "Level: ${level}${suffix}"
        metadata_tsv="../data/metadata_merged_l${level}${suffix}.tsv"
        taxonomy_tsv="../data/counts_merged_l${level}${suffix}.tsv"

        taxonomy_biom="../data/counts_merged_l${level}${suffix}.biom"
        conda deactivate && conda activate biom
        biom convert -i $taxonomy_tsv \
            -o $taxonomy_biom \
            --to-hdf5
        conda deactivate && conda activate qiime2-amplicon-2024.10

        taxonomy_qza="../data/counts_merged_l${level}${suffix}.qza"
        srun qiime tools import --type 'FeatureTable[Frequency]' \
             --input-path $taxonomy_biom \
             --output-path $taxonomy_qza

        bc_distances="../data/bc_distances_l${level}${suffix}.qza"
        srun qiime diversity beta --i-table $taxonomy_qza \
             --p-metric 'braycurtis' \
             --p-n-jobs 64 \
             --use-cache $cache_dir \
             --o-distance-matrix $bc_distances

        variables=("StudyID" "Primers" "HostSpecific" "Host" "Treatment" "RootCompartment")
        for variable in "${variables[@]}"; do
            echo "Running PERMANOVA for ${variable}…"
            srun qiime diversity adonis \
                 --i-distance-matrix $bc_distances \
                 --m-metadata-file $metadata_tsv \
                 --p-formula $variable \
                 --p-n-jobs 64 \
                 --use-cache $cache_dir \
                 --o-visualization "../data/permanova/${variable}_l${level}${suffix}.qzv"
            rm -rf "../data/permanova/${variable}_l${level}${suffix}"
            srun qiime tools extract \
                 --input-path "../data/permanova/${variable}_l${level}${suffix}.qzv" \
                 --output-path "../data/permanova/${variable}_l${level}${suffix}"
        done

        combined=("Primers+StudyID" "HostSpecific+StudyID" "Host+StudyID" "HostSpecific+RootCompartment+Treatment+StudyID" "HostSpecific+RootCompartment+Treatment+Primers+StudyID")
        for formula in "${combined[@]}"; do
            echo "Running PERMANOVA for ${formula}…"
            srun qiime diversity adonis \
                 --i-distance-matrix "$bc_distances" \
                 --m-metadata-file "$metadata_tsv" \
                 --p-formula $formula \
                 --p-n-jobs 64 \
                 --use-cache $cache_dir \
                 --o-visualization "../data/permanova/${formula}_l${level}${suffix}.qzv"
            rm -rf "../data/permanova/${formula}_l${level}${suffix}"
            srun qiime tools extract \
                 --input-path "../data/permanova/${formula}_l${level}${suffix}.qzv" \
                 --output-path "../data/permanova/${formula}_l${level}${suffix}"
        done
    done
done

