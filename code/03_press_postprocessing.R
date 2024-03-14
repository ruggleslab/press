###########################################################################
#
#                            prep_datasets
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: prep_datasets

# Load in the json from argparse
parser <- argparse::ArgumentParser()
parser$add_argument('--json', help = 'path to the config file')
args <- parser$parse_args()
params <- jsonlite::fromJSON(args$json)

# Output directory:
experiment <- "datasets"
outdir <- file.path(params$outdir, paste0(experiment))
dir.create(outdir, showWarnings = F)

# Would be nice to have one script that is going to test variious normalization and transformation methods

#======================== LIBRARIES ========================#
library(tidyverse)
library(readxl)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')



#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
