# This script is to generate json configs for the reduce press workflow

import os
import json

# Set up directories
outdir = "config"

#' The JSON config object contains three key-value pairs:
#' * "normalization": This key corresponds to the method used for normalization. In this case, "mor" is used.
#' * "geneset": This key corresponds to the path of the selected features file. In this case, the path is "output/feature_selection/rfe/logreg/selected_features.csv".
#' * "outdir": This key corresponds to the output directory where the results will be stored. In this case, the output directory is "output/reduce_press/rfeLOGREG_mor/".

# let's make a list of the normalization methods
normalization = ["mor", "tmm", "vst"]
geneset = []