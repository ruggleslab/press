#!/bin/bash
#SBATCH --job-name=press
#SBATCH --mem=32GB
#SBATCH --nodes=4
#SBATCH --array=1-1
#SBATCH -p cpu_short
#SBATCH --time=0-48:00:00
#SBATCH --output=logs/press_%A_%a.log

module load r/4.1.2

source /gpfs/data/ruggleslab/mm12865/mm12865_miniconda/etc/profile.d/conda.sh
conda activate main

# Get the parameters
configs=$(ls config/*.json)
json=$(echo $configs | cut -d " " -f ${SLURM_ARRAY_TASK_ID})

# Print the parameters
echo "Parameters"
echo "---------------------------------"
echo "json: ${json}"

echo " "
echo " "
echo " "

# copy the json file to the output directory

# Run prep_datasets.sh
echo "Running prep_datasets.R..."
echo "---------------------------------"
Rscript workflows/press/scripts/01_prep_datasets.R --json ${json}

echo " "
echo " "
echo " "

# Run press_model.sh
echo "Running press_model.py..."
echo "---------------------------------"
python workflows/press/scripts/02_press_model.py --json ${json}

echo " "
echo " "
echo " "

# # Run press postprocessing.R
echo "Running press_postprocessing.R..."
echo "---------------------------------"
Rscript workflows/press/scripts/03_press_postprocessing.R --json ${json}

echo " "
echo " "
echo " "

echo "Done!"