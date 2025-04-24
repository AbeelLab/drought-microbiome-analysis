#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
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

for level in $(seq 6 -1 2);
do
    echo "Level: ${level}"
    metadata_tsv="../data/merged_metadata_l${level}.tsv"
    taxonomy_tsv="../data/merged_only_taxa_l${level}.tsv"

    # Convert taxonomy to biom
    taxonomy_biom="../data/merged_only_taxa_l${level}.biom"
    conda deactivate && conda activate biom
    biom convert -i $taxonomy_tsv \
	 -o $taxonomy_biom \
	 --to-hdf5
    conda deactivate && conda activate qiime2-amplicon-2024.10

    taxonomy_qza="../data/merged_only_taxa_l${level}.qza"
    
    # Import taxonomic features as frequency table
    srun qiime tools import --type 'FeatureTable[RelativeFrequency]' \
    	 --input-path $taxonomy_biom \
	 --output-path $taxonomy_qza

    bc_distances="../data/bc_distances_l${level}.qza"
    # Create distance matrix
    srun qiime diversity beta --i-table $taxonomy_qza \
    	 --p-metric 'braycurtis' \
	 --p-n-jobs 32 \
	 --use-cache $cache_dir \
	 --o-distance-matrix $bc_distances

    # Run PERMANOVA
    # First run separately per variable
    # variables=("Study" "Primers" "Regions" "Treatment" "Host" "Location")
    # for variable in "${variables[@]}"; do
    # 	echo "Running PERMANOVA for ${variable}…"
    # 	srun qiime diversity adonis \
    #          --i-distance-matrix $bc_distances \
    #          --m-metadata-file $metadata_tsv \
    #          --p-formula $variable \
    #          --p-n-jobs 32 \
    # 	     --use-cache $cache_dir \
    #          --o-visualization "../data/permanova/${variable}_l${level}.qzv"
    # done

    # Then run sequential PERMANOVA for variables of interest
    # https://forum.qiime2.org/t/how-to-properly-design-formulas-in-adonis/17996/2
    # adonis performs a sequential test of terms which means the order in which you include your variables are actually important. Essentially, whatever variance is not explained by the first term is passed on to the second term, and then the remainder from that to the 3rd and so on
    combined=("Primers+Study" "Study+Primers")
    #("Host+Study" "Host+Primers" "Host+Regions" "Location+Study" "Location+Primers" "Location+Regions" "Host+Location+Study" "Host+Location+Regions" "Host+Location+Primers")
    for formula in "${combined[@]}"; do
    	echo "Running PERMANOVA for ${formula}…"
    	srun qiime diversity adonis \
             --i-distance-matrix "$bc_distances" \
             --m-metadata-file "$metadata_tsv" \
             --p-formula $formula \
             --p-n-jobs 32 \
	     --use-cache $cache_dir \
             --o-visualization "../data/permanova/${formula}_l${level}.qzv"
    done

    # Run differential abundance analysis
    
done
