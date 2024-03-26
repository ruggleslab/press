# This script is to generate json configs for the reduce press workflow

import os
import json
from glob import glob

# Set up directories
config_dir = "config"

#' The JSON config object contains three key-value pairs:
#' * "normalization": This key corresponds to the method used for normalization. In this case, "mor" is used.
#' * "geneset": This key corresponds to the path of the selected features file. In this case, the path is "output/feature_selection/rfe/logreg/selected_features.csv".
#' * "outdir": This key corresponds to the output directory where the results will be stored. In this case, the output directory is "output/reduce_press/rfeLOGREG_mor/".
#' * "model": This key corresponds to the model used for the PRESS algorithm. In this case, the model is "PRESS451".

# let's make a list of the normalization methods
normalization = ["mor", "tmm", "vst"]
geneset = glob("output/feature_selection/**/selected_features.csv", recursive=True)
model = ["PRESS451", "RF", "ADA", "LOGREG", "RUS"]

# make pairwise combinations of normalization and geneset
configs = []
for norm in normalization:
    for gene in geneset:
        for m in model:
            outdir = "output/reduce_press/" + gene.split("/")[2] + "_" + norm + "_" + m + "/"
            config = {
                "normalization": norm,
                "geneset": gene,
                "outdir": outdir,
                "model": m
            }
            configs.append(config)

# save each into the config directory
if not os.path.exists(outdir):
    os.makedirs(outdir)
for i, config in enumerate(configs):
    with open(f"{config_dir}/config_{i}.json", "w") as f:
        json.dump(config, f, indent=4)