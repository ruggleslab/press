#!/bin/bash

# if the log exists wipe it
# if [ -f logs/press_model_$(date +"%Y%m%d").log ]; then
#     rm logs/press_model_$(date +"%Y%m%d").log
# fi
exec 1>>logs/press_model_$(date +"%Y%m%d").log 2>&1

# Get the parameters
json=$1

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
Rscript workflow/steps/press/01_prep_datasets.R --json ${json}

echo " "
echo " "
echo " "

# Run press_model.sh
echo "Running press_model.py..."
echo "---------------------------------"
python workflow/steps/press/02_press_model.py --json ${json}

echo " "
echo " "
echo " "

# # Run press postprocessing.R
# echo "Running press_postprocessing.R..."
# echo "---------------------------------"
# Rscript workflow/steps/press/03_press_postprocessing.R --json ${json}


# echo "Script execution completed."
