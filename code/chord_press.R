###########################################################################
#
#                            chord_press
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-09-20
# Script Name: chord_press
# Output directory:
experiment <- "chord_press"
outdir <- file.path("output", experiment)
dir.create(outdir, showWarnings = FALSE)

#======================== LIBRARIES ========================
library(tidyverse)
library(glue)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/signature_functions.R")

#======================== CODE ========================
meta1 <- read.csv("/Users/muller/Library/CloudStorage/GoogleDrive-mm12865@nyu.edu/My Drive/RugglesLab/datasets/chord/processed/chord_press_scores.csv", row.names = 1)
meta2 <- read.csv("/Users/muller/Library/CloudStorage/GoogleDrive-mm12865@nyu.edu/My Drive/RugglesLab/datasets/chord/preprocessing/metadata.csv", row.names = 1)
meta <- merge(meta1, meta2, by = "row.names")

#======================== Comparisons ========================
# get the variables for a1c, ldl, and hscrp
vars <- c("bl_a1c", "bl_ldl", "bl_crp")

res <- meta %>%
    filter(Study_Timepoint == 1) %>%
    compare_one_to_many("scores", vars, file.path(outdir, "chord_press"), plot = TRUE)


















#======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))