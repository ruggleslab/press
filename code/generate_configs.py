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

# let's make a list of the normalization methods
normalization = ["mor", "tmm", "vst"]
geneset = glob("output/feature_selection/**/selected_features.csv", recursive=True)

# make pairwise combinations of normalization and geneset
configs = []
for norm in normalization:
    for gene in geneset:
        outdir = "output/reduce_press/" + gene.split("/")[2] + "_" + norm + "/"
        config = {
            "normalization": norm,
            "geneset": gene,
            "outdir": outdir
        }
        configs.append(config)

# save each into the config directory
if not os.path.exists(outdir):
    os.makedirs(outdir)
for i, config in enumerate(configs):
    with open(f"{config_dir}/config_{i}.json", "w") as f:
        json.dump(config, f, indent=4)
        print(f"config_{i}.json saved into {config_dir}")