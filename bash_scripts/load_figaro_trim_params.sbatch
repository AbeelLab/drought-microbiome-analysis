#!/bin/sh

#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --qos=short
#SBATCH --time=00:05:00
#SBATCH --mem=4GB
#SBATCH --mail-type=END
#SBATCH --output=slurm_outputs/%x_%j.out
#SBATCH --error=slurm_outputs/%x_%j.out

/usr/bin/scontrol show job -d "$SLURM_JOB_ID"

data_path=$1
study_id=$2
num_batches=$3

source /tudelft.net/staff-umbrella/abeellab/bmcosma/miniconda3/etc/profile.d/conda.sh

study_dir="${data_path}/${study_id}"
working_dir="${study_dir}/data_batches"

figaro_ran_succesfully=true
i=1
while [ $i -le $num_batches ]
do
    qiime_dir="${working_dir}/qiime-dir-batch${i}/figaro-output"
    if [ ! -f "${qiime_dir}/trimParameters.json" ]; then
        figaro_ran_successfully=false
        break
    fi
    i=$(( i + 1 ))
done

if [ "$figaro_ran_successfully" = true ]; then
    conda activate general
    cd ../python_scripts/
    python3 extract_figaro_parameters.py --study $study_id
    cd ../bash_scripts/
    conda deactivate
else
    echo "FIGARO did not run successfully; will truncate at 90% amplicon length"
fi
