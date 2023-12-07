###########################################################################
#
#                            scratch
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-07-28
# Script Name: scratch
# Output directory:
experiment <- "scratch"
run <- 2
outdir <- file.path('output', paste0(experiment, '__run_', run))
dir.create(outdir, showWarnings = F)

# Notes about the experiment run:
notes <- "scratch book" #nolint

# save notes to file
write(notes, file.path(outdir, "notes.txt"))

#======================== LIBRARIES ========================#
packages <- c("lintr", "httpgd", "languageserver", "devtools", "sys", "dplyr", "tidyverse", "ggplot2", "ggpubr", "SummarizedExperiment")
pkgs <- lapply(packages, function(x) suppressMessages(require(x, character.only=T,quietly=T))) # nolint
print(packages[!sapply(pkgs, isTRUE)])

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R')

#======================== CODE ========================#

scores <- read.csv('export.csv', header = T, stringsAsFactors = F)
scores <- scores[-84, ]
scores$scores <- scale(scores$preds)
hist(scores$scores)

# save it
write.csv(scores, file.path(outdir, "scores.csv"), row.names = F)















#======================== END ========================#
save.image(file.path(outdir, "image.RData"))

sink(file.path(outdir, "session.log"), append = TRUE, split = TRUE)
sessionInfo()
sink()
