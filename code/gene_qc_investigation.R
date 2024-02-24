###########################################################################
#
#                            gene_qc_investigation
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-01-21
# Script Name: gene_qc_investigation
# Output directory:
experiment <- "gene_qc_investigation"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

# Notes about the experiment run:
notes <- "So let's choose some genes of interest and look at those since shit is kinda not lining up for chord..." #nolint

# save notes to file
write(notes, file.path(outdir, "notes.txt"))

#======================== LIBRARIES ========================#
library(lintr) #nolint
library(httpgd) #nolint
library(languageserver) #nolint
library(devtools)
library(tidyverse)
library(glue)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')

#======================== DATA ========================#
# load in the press genes
pace_dge <- read.csv('data/PACE_hyper_v_hypo_deseqoutput.csv', header = TRUE, row.names = 1)
duke_dge <- read.csv('data/duke_validation_run3/dge_analysis/comp_group1__hyper_v_nothyper_AGCONTROL_deseqout.csv', header = TRUE, row.names = 1)

# comparisons of interest
press <- read.csv('data/clean/press_genes.csv', header = FALSE) %>% pull(.)
high_low_comps <- c('HHvLL', 'HH.HMvML.LL', 'Universal_HvL', 'UHvLL', 'HHvUL')

##### CHORD #####
chord_counts <- read.csv('output/chord_data_cleaning/chord_press_counts.csv', row.names = 1)
chord_scores <- read.csv('output/chord_press_scores/chord_data_with_press_scores.csv', row.names = 1, na.strings = 'nan')

##### DUKE #####
duke_counts <- read.csv('output/duke_all_data/duke_all_data.csv', row.names = 1)
duke_scores <- read.csv('output/duke_all_data/duke_all_data_scores.csv', row.names = 1)

##### PACE #####
pace_counts <- read.csv('output/PACE_GEO_SUBM/press_counts.csv', row.names = 1)
pace_scores <- read.csv('output/pace_press_scoring_tiles__run_4/pace_press_metadata_with_scores_and_tiles.csv')

#### HARP ####
# TODO: add harp data








#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
