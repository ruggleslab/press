###########################################################################
#
#                            pace_scoring_investigation
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-06-09
# Script Name: pace_scoring_investigation
# Notes:
# This script looks into the scoring we have for the pace project.

# Output directory:
experiment <- 'pace_scoring_investigation'
run <- 1
outdir <- file.path('output', paste0(experiment, '__run_', run))
dir.create(outdir, showWarnings = F)

#======================== LIBRARIES ========================#
packages <- c("devtools", "sys", "dplyr", "tidyverse", "ggplot2", "ggpubr", "SummarizedExperiment", "readxl", "cowplot")
for (pkg in packages) {library(pkg, character.only = T, quietly = T)}

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R')

#======================== CODE ========================#
# load in the pace_scoring data
pace_scoring <- read_excel('data/pace_scoring.xlsx', na = 'NA')

pace_scoring <- pace_scoring %>%
    mutate(labels = as.factor(labels))

# get the median, 5th percentile, 25th percentile, 75th percentile, and 95th percentile
pace_scoring_summary <- pace_scoring %>%
    summarise(
        median = median(scores),
        p5 = quantile(scores, probs = 0.05),
        p25 = quantile(scores, probs = 0.25),
        p75 = quantile(scores, probs = 0.75),
        p95 = quantile(scores, probs = 0.95)
    )

# save the summary data
write.csv(pace_scoring_summary, file.path(outdir, 'pace_cohort_scoring_summary.csv'), row.names = F)

# make a histogram of the data
score_histogram <- ggplot(pace_scoring, aes(x = scores)) +
    geom_histogram(bins = 50, fill = 'blue', alpha = 0.75) +
    labs(title = 'Histogram of Pace Scoring (All Samples)', x = 'Score', y = 'Count') +
    theme_classic2() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        # increase legend text size
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)
    ) +
    # add density plot
    geom_density(
        aes(y = ..density..),
        fill = 'blue'
    ) +
    # add in the median, 5th percentile, 25th percentile, 75th percentile, and 95th percentile
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = median, color = 'median'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p5, color = 'p5'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p25, color = 'p25'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p75, color = 'p75'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p95, color = 'p95'),
        linetype = 'dashed',
        size = 1
    )

# get the median, 5th percentile, 25th percentile, 75th percentile, and 95th percentile for the derivation data
pace_scoring_summary <- pace_scoring %>%
    filter(!is.na(labels)) %>%
    summarise(
        median = median(scores),
        p5 = quantile(scores, probs = 0.05),
        p25 = quantile(scores, probs = 0.25),
        p75 = quantile(scores, probs = 0.75),
        p95 = quantile(scores, probs = 0.95)
    )

# save the summary data
write.csv(pace_scoring_summary, file.path(outdir, 'pace_derivation_cohort_scoring_summary.csv'), row.names = F)

# make a histogram of the derivation data
score_derivation_histogram <- ggplot(filter(pace_scoring, !is.na(labels)), aes(x = scores)) +
    geom_histogram(bins = 50, fill = 'blue', alpha = 0.75) +
    labs(title = 'Histogram of Pace Scoring (Derivation Cohort)', x = 'Score', y = 'Count') +
    theme_classic2() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        # increase legend text size
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 24)
    ) +
    # add density plot
    geom_density(
        aes(y = ..density..),
        fill = 'blue',
        alpha = 0.75,
    ) +
    # add in the median, 5th percentile, 25th percentile, 75th percentile, and 95th percentile
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = median, color = 'median'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p5, color = 'p5'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p25, color = 'p25'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p75, color = 'p75'),
        linetype = 'dashed',
        size = 1
    ) +
    geom_vline(
        data = pace_scoring_summary,
        aes(xintercept = p95, color = 'p95'),
        linetype = 'dashed',
        size = 1
    )

# plot the histograms side by side
pace_histograms <- plot_grid(score_histogram, score_derivation_histogram, nrow = 2)

# save the histograms
ggsave(file.path(outdir, 'pace_histograms.png'), pace_histograms, width = 20, height = 16, units = 'in', dpi = 300)

# make a boxplot of the data
score_boxplot <- ggplot(pace_scoring, aes(x = labels, y = scores, fill = labels)) +
    geom_boxplot() +
    labs(title = 'Boxplot of Pace Scoring', x = 'Labels (0:Normal, 1:Hyper)', y = 'Score') +
    theme_classic2() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    )

# save the boxplot
ggsave(file.path(outdir, 'pace_boxplot.png'), score_boxplot, width = 20, height = 16, units = 'in', dpi = 300)


#======================== SINGSCORES PACE VERSUS DUKE ========================#
library(singscore)

# load in the up and down gene sets
down_genes <- read.table('data/custom_mgc_hyper_down.txt', header = F, stringsAsFactors = F) %>% pull()
up_genes <- read.table('data/custom_mgc_hyper_up.txt', header = F, stringsAsFactors = F) %>% pull()

# load in the pace data
pace_counts <- read.csv('data/hyper_feature_outtable.csv', row.names = 1)
# this is the useful one!!
pace_counts_raw <- read.table('data/plt_filtrawcounttab.txt', sep = '\t', row.names = 1, header = T)

# subset the raw counts to only have samples in pace_counts
pace_counts_raw <- pace_counts_raw[, rownames(pace_counts)]
# remove the _ from the sample names
colnames(pace_counts_raw) <- gsub('_', '', colnames(pace_counts_raw))
rownames(pace_counts) <- gsub('_', '', colnames(pace_counts_raw))


# load in the duke data
load('data/summarized_experiments/se_duke_pltRNA.rdata')
se_duke_plt

# get the duke counts
duke_counts <- assay(se_duke_plt)

# duke_counts <- duke_counts %>% add_missing_rows(c(down_genes, up_genes))

# get the singscores for the pace data
pace_singscores <- simpleScore(rankGenes(pace_counts_raw), up_genes, down_genes)

# add the singscores to the pace_scoring data
pace_scoring_derivation_cohort <- pace_scoring %>%
    filter(!is.na(labels))

# add the singscores to the pace_scoring data
pace_scoring_derivation_cohort <- pace_scoring_derivation_cohort %>%
    mutate(singscore = pace_singscores)

# get the singscores for the duke data
duke_singscores <- simpleScore(rankGenes(duke_counts), up_genes, down_genes)

# add the singscores to the duke data
duke_scoring <- se_duke_plt %>% 
    colData() %>%
    as.data.frame() %>%
    rownames_to_column('sample') %>%
    mutate(singscore = duke_singscores)

# plot the singscores
singscore_plot_pace <- ggplot(pace_scoring_derivation_cohort, aes(x = labels, y=singscore$TotalScore, fill = labels)) +
    geom_boxplot() +
    labs(title = 'Boxplot of Pace Singscores', x = 'Labels (0:Normal, 1:Hyper)', y = 'Singscore') +
    theme_classic2() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    ) +
    stat_compare_means(label.y = 0.5, size = 6)

singscore_plot_duke <- ggplot(duke_scoring, aes(x = Phenotype, y=singscore$TotalScore, fill = Phenotype)) +
    geom_boxplot() +
    labs(title = 'Boxplot of Duke Singscores', x = 'Labels (0:Normal, 1:Hyper)', y = 'Singscore') +
    theme_classic2() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    ) +
    stat_compare_means(label.y = 0.5, size = 6)

# plot the singscores side by side
singscore_plots <- plot_grid(singscore_plot_pace, singscore_plot_duke, nrow = 2)

# save the singscores
ggsave(file.path(outdir, 'singscores.png'), singscore_plots, width = 12, height = 14, units = 'in', dpi = 300)


# correlation of singscores to scores
cor_plot_pace <- ggplot(pace_scoring_derivation_cohort, aes(x = scores, y=singscore$TotalScore)) +
    geom_point() + geom_smooth(method = 'lm') +
    labs(title = 'Correlation of Singscores to Model Scores in Pace', x = 'Model Scores', y = 'Singscore Scores') +
    theme_classic2() +
    theme(
        plot.title = element_text(hjust = 0.5, size = 24),
        axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)
    ) +
    stat_cor(method = 'spearman', label.x = 0.5, label.y = 0.5, size = 12)

ggsave(file.path(outdir, 'pace_singscore_correlation.png'), cor_plot_pace, width = 20, height = 16, units = 'in', dpi = 300)