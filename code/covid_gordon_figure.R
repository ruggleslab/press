###########################################################################
#
#                            covid_gordon_figure
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-07-29
# Script Name: covid_gordon_figure
# Output directory:
experiment <- "covid_gordon_figure"
outdir <- file.path("output", experiment)
dir.create(outdir, showWarnings = FALSE)

#======================== LIBRARIES ========================
library(tidyverse)
library(glue)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/enrichment_functions.R")

#======================== CODE ========================
# So we are going to just make a figure showing that PRESS works
# in COVID for JBs presentation at the Gordon conference
# this can be done with the data we have already
# we will want to make a 3 panels --
# 1. the PRESS score for our COVID data
# 2. the PRESS enrichment GSEA for our COVID data
# 3. the PRESS enrichment GSEA for the COVID Utah Data

press <- pull(read.csv('data/press451_genes.csv', header = T))

#======================== Panel 1 ========================
dir.create(file.path(outdir, "covid_scores"), showWarnings = FALSE)
# press for our COVID data
counts_raw <- read.csv("data/20240729_covid_figure/covid_covidcontrol/intermediate_files/covid_covidcontrol__rawcounttable.csv", row.names = 1)
meta <- read.csv("data/20240729_covid_figure/covid_covidcontrol/intermediate_files/covid_covidcontrol__metatable.csv", row.names = 1)
rownames(meta) <- gsub("-", "\\.", rownames(meta))

missing_genes <- setdiff(press, rownames(counts_raw))
counts_raw_clean <- counts_raw %>% 
    add_missing_rows(press) %>%
    filter(rownames(.) %in% press)

se <- make_se(counts_raw_clean, meta)
dds <- DESeqDataSet(se, design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_logmor <- t(normalize_counts(dds, "mor", log = TRUE))
write.csv(counts_logmor, file.path(outdir, "covid_scores", "press_counts.csv"))

dat <- glue("{outdir}/covid_scores/press_counts.csv")
out <- glue("{outdir}/covid_scores/press_scores.csv")
system(glue("python3 code/run_press.py --data {dat} --out {out}"))

scores <- read.csv(file.path(outdir, "covid_scores", "press_scores.csv"), row.names = 1)
dat <- cbind(meta, scores)

covid_press_plot <- ggplot(dat, aes(x = factor(comp_covid__covid_v_control), y = scores, fill = factor(comp_covid__covid_v_control))) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "COVID vs Control", y = "PRESS Score") +
    geom_hline(yintercept = 0.38, linetype = "dashed") +
    annotate("text", x = 0.8, y = 0.38, label = "Hyperreactive", vjust = -1, color = "red", size = 6) +
    annotate("text", x = 0.8, y = 0.38, label = "Normoreactive", vjust = 1, color = "blue", size = 6) +
    theme(legend.position = "none")
ggsave(file.path(outdir, "covid_scores", "covid_press_plot.pdf"), covid_press_plot, width = 5, height = 5)











#======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
