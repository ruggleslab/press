#!/bin/bash
#SBATCH --job-name=press
#SBATCH --mem=16GB
#SBATCH --nodes=1
#SBATCH -p cpu_short
#SBATCH --time=0-48:00:00
#SBATCH --output=logs/reduce_press_%j.log


# WORKING DIRECTORY
parentdir="$(dirname "dir")" # get the parent directory
cd parentdir


# SCRIPTS

# end of script

##!/bin/bash

# This script requires the use of:
#   - config/*.json
#   - press_model.sh
#       - 01_datasets.R
#       - 02_press_model.R
#       - 03_press_postprocessing.R

#' The JSON config object contains three key-value pairs:
#' * "normalization": This key corresponds to the method used for normalization. In this case, "mor" is used.
#' * "geneset": This key corresponds to the path of the selected features file. In this case, the path is "output/feature_selection/rfe/logreg/selected_features.csv".
#' * "outdir": This key corresponds to the output directory where the results will be stored. In this case, the output directory is "output/reduce_press/rfeLOGREG_mor/".

# make a list of the configs in the config directory

echo "Number of configs: $(echo $configs | wc -w)"
echo " "
echo " "

# loop through the configs and run the press_model.sh script
for config in $configs; do
    echo "Running reduce_press.sh..."
    echo "---------------------------------"
    echo "config: ${config}"
    echo " "
    echo Sbatching press_model.sh ${config}
    sbatch workflows/steps/press/press_model.sh ${config}
    echo " "
    echo " "
    echo " "
done