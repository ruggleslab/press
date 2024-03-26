###########################################################################
#
#                            press_therapy_effects
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-03-20
# Script Name: press_therapy_effects
# Output directory:
experiment <- "press_therapy_effects"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

#======================== LIBRARIES ========================#
library(tidyverse)
library(glue)
library(rstatix)
library(ggpubr)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')

#======================== CODE ========================#
# So as a follow up to PRESS it would be ideal if we can look at what sort of therapies have an effect on press. We can do this in:
# - PACE
# - DUKE
# - Treated MKs
# - Treated PLTs

# lets start by loading in the data
treated_dds_list <- readRDS('data/treated_cell_lines_data/dataset_list.rds')
datasets <- c("lupus", "megakaryocyte", "agonist", "ifn")
press451 <- pull(read.csv('data/press451_genes.csv'))

dir.create(file.path(outdir, "datasets"), showWarnings = F)
dir.create(file.path(outdir, "press_scores"), showWarnings = F)
score_list <- list()
for (dat in datasets) {
    se <- treated_dds_list[[dat]]
    data_norm <- assay(se) %>% add_missing_rows(press451)
    # dds <- DESeqDataSet(se, ~1)
    # data_norm <- normalize_counts(dds, method = "log2-mor")
    data_norm <- t(data_norm[press451,])
    data <- glue('{outdir}/datasets/{dat}.csv')
    write.csv(data_norm, data)
    out <- glue('{outdir}/press_scores/{dat}_press.csv')
    system(glue('python3 code/run_press.py --data {data} --out {out}'))

    scores <- read.csv(out, row.names = 1)
    meta <- colData(se)
    meta_with_scores <- cbind(meta, scores)
    score_list[[dat]] <- meta_with_scores
}

#======================== Dataset Subanalyses ========================
## Lupus
dir.create(file.path(outdir, "lupus"), showWarnings = F)
lupus <- as.data.frame(score_list$lupus)
write.csv(lupus, file.path(outdir, "lupus", "lupus_meta.csv"))
# TODO

## MKs
dir.create(file.path(outdir, "megakaryocyte"), showWarnings = F)
mks <- as.data.frame(score_list$megakaryocyte)
mks$antiplt <- factor(
    mks$antiplt, 
    levels = c("PBS", "Aspirin", "DMSO", "AZD1283", "Aspirin_AZD1283"),
    labels = c("PBS", "ASA", "DMSO", "P2Y12", "ASA+P2Y12")
    )
write.csv(mks, file.path(outdir, "megakaryocyte", "mks_meta.csv"))
hist <- ggplot(mks, aes(x = scores)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0.38, color = "red") +
    theme_bw() +
    labs(title = "Megakaryocyte Press Scores", x = "Press Score", y = "Counts")
ggsave(file.path(outdir, "megakaryocyte", "mks_hist.png"), hist, width = 5, height = 5)
boxplot <- ggplot(mks, aes(x = antiplt, y = scores)) +
    geom_boxplot() +
    geom_point(aes(color = donor)) +
    geom_line(aes(group = donor, color = donor)) +
    geom_hline(yintercept = 0.38, color = "red", linetype = "dashed") +
    theme_bw() +
    labs(title = "Megakaryocyte Press Scores", x = NULL, y = "Press Score")
ggsave(file.path(outdir, "megakaryocyte", "mks_boxplot.png"), boxplot, width = 5, height = 5)
head(mks)
unpaired_stats <- mks %>% t_test(scores ~ antiplt, paired = F)
write.csv(unpaired_stats, file.path(outdir, "megakaryocyte", "mks_unpaired_stats.csv"))
paired_stats <- mks %>% 
    filter(antiplt %in% c("DMSO", "P2Y12", "ASA+P2Y12")) %>%
    arrange(donor, antiplt) %>%
    mutate(antiplt = factor(antiplt, levels = c("DMSO", "P2Y12", "ASA+P2Y12"))) %>%
    t_test(scores ~ antiplt, paired = T, detailed = T)
write.csv(paired_stats, file.path(outdir, "megakaryocyte", "mks_paired_stats.csv"))


## Agonist
dir.create(file.path(outdir, "agonist"), showWarnings = F)
agonist <- as.data.frame(score_list$agonist)
agonist$agonists <- factor(agonist$agonists, levels = c("PBS", "Thrombin", "ADP", "Rhodocytin", "Collagen", "DMSO", "U46619"))
hist <- ggplot(agonist, aes(x = scores)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0.38, color = "red") +
    theme_bw() +
    labs(title = "Agonist Press Scores", x = "Press Score", y = "Counts")
ggsave(file.path(outdir, "agonist", "agonist_hist.png"), hist, width = 5, height = 5)
comps <- list(
    c('PBS', 'Collagen'),
    c('PBS', 'Rhodocytin'),
    c('PBS', 'ADP'),
    c('PBS', 'Thrombin'),
    c('DMSO', 'U46619')
)
boxplot <- ggplot(agonist, aes(x = agonists, y = scores)) +
    geom_boxplot() +
    geom_point(aes(color = donor)) +
    geom_line(aes(group = donor, color = donor)) +
    geom_hline(yintercept = 0.38, color = "red", linetype = "dashed") +
    stat_compare_means(method = 't.test', comparisons = comps, paired = T) +
    theme_bw() +
    labs(title = "Agonist Press Scores", x = NULL, y = "Press Score", subtitle = "Paired T-Test")
ggsave(file.path(outdir, "agonist", "agonist_boxplot.png"), boxplot, width = 5, height = 5)
unpaired_stats <- agonist %>% t_test(scores ~ agonists)
write.csv(unpaired_stats, file.path(outdir, "agonist", "agonist_unpaired_stats.csv"))

## IFN
dir.create(file.path(outdir, "ifn"), showWarnings = F)
ifn <- as.data.frame(score_list$ifn)
ifn$IFNs <- factor(ifn$IFNs, levels = c("PBS", "IFNalpha", "IFNgamma", "LPS", "Poly"))
write.csv(ifn, file.path(outdir, "ifn", "ifn_meta.csv"))
hist <- ggplot(ifn, aes(x = scores)) +
    geom_histogram(bins = 30) +
    geom_vline(xintercept = 0.38, color = "red") +
    theme_bw() +
    labs(title = "IFN Press Scores", x = "Press Score", y = "Counts")
ggsave(file.path(outdir, "ifn", "ifn_hist.png"), hist, width = 5, height = 5)
comps <- list(
    c('PBS', 'IFNalpha'),
    c('PBS', 'IFNgamma'),
    c('PBS', 'LPS'),
    c('PBS', 'Poly')
)
boxplot <- ggplot(ifn, aes(x = IFNs, y = scores)) +
    geom_boxplot() +
    geom_point(aes(color = donor)) +
    geom_line(aes(group = donor, color = donor)) +
    geom_hline(yintercept = 0.38, color = "red", linetype = "dashed") +
    stat_compare_means(method = 't.test', comparisons = comps, paired = T) +
    theme_bw() +
    labs(title = "IFN Press Scores", x = NULL, y = "Press Score", subtitle = "Paired T-Test")
ggsave(file.path(outdir, "ifn", "ifn_boxplot.png"), boxplot, width = 5, height = 5)




#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
