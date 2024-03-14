##!/bin/bash

# This script requires the use of:
#   - config/*.json
#   - press_model.sh
#       - 01_datasets.R
#       - 02_press_model.R
#       - 03_press_postprocessing.R

# if the log exists wipe it
if [ -f logs/reduce_press_$(date +"%Y%m%d").log ]; then
    rm logs/reduce_press_$(date +"%Y%m%d").log
fi
exec 1>>logs/reduce_press_$(date +"%Y%m%d").log 2>&1

# make a list of the configs in the config directory
configs=$(ls config/*.json)
echo "Number of configs: $(echo $configs | wc -w)"
echo " "
echo " "

# loop through the configs and run the press_model.sh script
for config in $configs; do
    echo "Running reduce_press.sh..."
    echo "---------------------------------"
    echo "config: ${config}"
    echo " "
    echo code/press_model.sh ${config}
    code/press_model.sh ${config}
    echo " "
    echo " "
    echo " "
done