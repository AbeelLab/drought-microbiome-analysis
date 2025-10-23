#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=short
#SBATCH --time=02:00:00
#SBATCH --mem=32GB
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

# Set up QIIME2 cache
# cache_dir="/tmp/job_${SLURM_JOB_ID}"
# srun qiime tools cache-create --cache $cache_dir
# export TMPDIR=$cache_dir

# # Clean up cache after job execution
# function clean_up {
#     echo "Before cleanup:"
#     find /tmp -user bmcosma -print0 | du -chs --files0-from=-
#     echo "---"
    
#     rm -rf $cache_dir

#     echo "After cleanup:"
#     find /tmp -user bmcosma -print0 | du -chs --files0-from=-
#     exit
# }

# trap 'clean_up' EXIT

for level in 6; do
    for bc in "" "_batch_corrected"; do
        for compartment in "bulk_soil" "endosphere" "rhizosphere"; do
	    echo "Level: ${level}${suffix}"
      	    metadata_tsv="../datasets/${compartment}/metadata.tsv"
	    taxonomy_tsv="../datasets/${compartment}/counts_genus_for_signature${bc}.tsv"
	    taxonomy_tsv_t="${taxonomy_tsv%.tsv}.transposed.tsv"

	    taxonomy_biom="../raw_data/${compartment}_counts_genus_for_signature${bc}.biom"
	    conda deactivate && conda activate biom
	    biom convert -i $taxonomy_tsv_t \
		 -o $taxonomy_biom \
		 --to-hdf5
	    conda deactivate && conda activate qiime2-amplicon-2024.10

	    taxonomy_qza="../raw_data/${compartment}_counts_genus_for_signature${bc}.qza"
	    srun qiime tools import --type 'FeatureTable[Frequency]' \
		 --input-path $taxonomy_biom \
		 --output-path $taxonomy_qza

	    bc_distances="../raw_data/${compartment}_counts_genus_for_signature${bc}_distances.qza"
	    srun qiime diversity beta --i-table $taxonomy_qza \
		 --p-metric 'braycurtis' \
		 --p-n-jobs 64 \
		 --o-distance-matrix $bc_distances

	    variables=("StudyID" "Primers" "HostSpecific" "Host" "Treatment")
	    for variable in "${variables[@]}"; do
		echo "Running PERMANOVA for ${variable}…"
		srun qiime diversity adonis \
		     --i-distance-matrix $bc_distances \
		     --m-metadata-file $metadata_tsv \
		     --p-formula $variable \
		     --p-n-jobs 16 \
		     --verbose \
		     --o-visualization "../raw_data/permanova/${variable}_l${level}_${compartment}${bc}.qzv"
		rm -rf "../raw_data/permanova/${variable}_l${level}_${compartment}${bc}"
		srun qiime tools extract \
		     --input-path "../raw_data/permanova/${variable}_l${level}_${compartment}${bc}.qzv" \
		     --output-path "../raw_data/permanova/${variable}_l${level}_${compartment}${bc}"
	    done
	done
    done
done

