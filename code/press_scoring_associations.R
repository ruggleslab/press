###########################################################################
#
#                            pace_press_scoring_tiles
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-06-12
# Script Name: press_scoring_associations
# Notes:
# Now that we have have a good model, we can check to see

# Output directory:
experiment <- "press_scoring_associations"
outdir <- file.path('output', experiment)
dir.create(outdir, showWarnings = F)

#======================== LIBRARIES ========================#
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(SummarizedExperiment)
library(readxl)
library(cowplot)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')

#======================== PACE ========================#
# load in the scores
# press_scores <- read_xlsx('docs/pace_scoring_prior_set.xlsx')
press_scores <- read_xlsx('data/pace_scoring.xlsx')

# set the row names to the subject ids
press_scores <- press_scores %>%
    as.data.frame() %>%
    column_to_rownames('Sample_ID')


# load in the metadata
metadata <- read.csv('data/pace/comb_updated2021.csv', stringsAsFactors = T, header = T)
metadata <- metadata %>% 
    column_to_rownames('name_updated')

# get the intersection of the two
matched_samples <- intersect(rownames(press_scores), rownames(metadata))

# add the scores to the metadata
metadata <- cbind(metadata[matched_samples,], press_scores[matched_samples,])

# add 3 and 4 tiles based on the scores to the metadata
metadata <- metadata %>% 
    mutate(
        # add tiles
        tile_2 = as.factor(ntile(scores, 2)),
        tile_3 = as.factor(ntile(scores, 3)),
        tile_4 = as.factor(ntile(scores, 4))
        ) %>%
        # make censor variables factors
        mutate_at(vars(starts_with('censor')), as.factor)

# save the metadata
write.csv(metadata, file.path(outdir, 'pace_press_metadata_with_scores_and_tiles.csv'))

# make histograms of the scores for each tile
# 2 tiles
p2 <- metadata %>% 
    ggplot(aes(x = scores, fill = tile_2)) +
    geom_histogram() +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    labs(x = 'PACE Press Score', y = 'Count', title = 'PACE Press Scores by 2 Tiles')
# 3 tiles
p3 <- metadata %>% 
    ggplot(aes(x = scores, fill = tile_3)) +
    geom_histogram() +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    labs(x = 'PACE Press Score', y = 'Count', title = 'PACE Press Scores by 3 Tiles')
# 4 tiles
p4 <- metadata %>% 
    ggplot(aes(x = scores, fill = tile_4)) +
    geom_histogram() +
    theme_bw() +
    theme(text = element_text(size = 20)) +
    labs(x = 'PACE Press Score', y = 'Count', title = 'PACE Press Scores by 4 Tiles')
# combine the plots
hist_plots <- plot_grid(p2, p3, p4, nrow = 3)
# save the plots
ggsave(file.path(outdir, 'pace_press_scoring_tiles_histograms.pdf'), hist_plots, width = 10, height = 10)


# make barplots of the scores for each tile based upon:
# - censor_MACLE2
# - censor_MALE2
# - censor_MACE

# censor_MACLE2
p2 <- metadata %>% 
    ggplot(aes(x = tile_2, fill = censor_MACLE2)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '2 Tiles', y = 'Count', title = 'PACE Press Scores by 2 Tiles and censor_MACLE2')
p3 <- metadata %>% 
    ggplot(aes(x = tile_3, fill = censor_MACLE2)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '3 Tiles', y = 'Count', title = 'PACE Press Scores by 3 Tiles and censor_MACLE2')
p4 <- metadata %>% 
    ggplot(aes(x = tile_4, fill = censor_MACLE2)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '4 Tiles', y = 'Count', title = 'PACE Press Scores by 4 Tiles and censor_MACLE2')
# combine the plots
macle_plots <- plot_grid(p2, p3, p4, nrow = 3)
# save the plots
ggsave(file.path(outdir, 'pace_press_scoring_tiles_censor_MACLE2.pdf'), macle_plots, width = 10, height = 10)

# censor_MALE2
p2 <- metadata %>% 
    ggplot(aes(x = tile_2, fill = censor_MALE2)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '2 Tiles', y = 'Count', title = 'PACE Press Scores by 2 Tiles and censor_MALE2')
p3 <- metadata %>% 
    ggplot(aes(x = tile_3, fill = censor_MALE2)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '3 Tiles', y = 'Count', title = 'PACE Press Scores by 3 Tiles and censor_MALE2')
p4 <- metadata %>%
    ggplot(aes(x = tile_4, fill = censor_MALE2)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '4 Tiles', y = 'Count', title = 'PACE Press Scores by 4 Tiles and censor_MALE2')
# combine the plots
male2_plots <- plot_grid(p2, p3, p4, nrow = 3)
# save the plots
ggsave(file.path(outdir, 'pace_press_scoring_tiles_censor_MALE2.pdf'), male2_plots, width = 10, height = 10)

# censor_MACE
p2 <- metadata %>% 
    ggplot(aes(x = tile_2, fill = censor_MACE)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '2 Tiles', y = 'Count', title = 'PACE Press Scores by 2 Tiles and censor_MACE')
p3 <- metadata %>% 
    ggplot(aes(x = tile_3, fill = censor_MACE)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '3 Tiles', y = 'Count', title = 'PACE Press Scores by 3 Tiles and censor_MACE')
p4 <- metadata %>%
    ggplot(aes(x = tile_4, fill = censor_MACE)) +
    geom_bar(position = 'dodge') +
    theme_matt(18) +
    labs(x = '4 Tiles', y = 'Count', title = 'PACE Press Scores by 4 Tiles and censor_MACE')
# combine the plots
mace_plots <- plot_grid(p2, p3, p4, nrow = 3)
# save the plots
ggsave(file.path(outdir, 'pace_press_scoring_tiles_censor_MACE.pdf'), mace_plots, width = 10, height = 10)

#========= Correlation Plots =========#
# get the correlation between the scores and the epi variables
scores <- metadata$scores
epi_variables <- metadata %>%
    dplyr::select(
        contains(c('180s', '300s', 'max')), 
        scores
        ) %>%
    dplyr::select(
        starts_with(
            c(
                'epi_20', 'epi_01', 'epi_04',
                'adp_20', 'adp_01', 'adp_04',
                'col_10', 'col_02',
                'ser10', 
                'aa_1600', 'aaexvivo_1600'
                )
            )
        )

### ALL SAMPLES
# map cor.test to get a p-value and correlation
epi_cor <- map_df(
    epi_variables, 
    ~cor.test(.x, scores) %>% 
        broom::tidy() %>% 
        dplyr::select(p.value, estimate) %>% 
        dplyr::rename(p = p.value, scores = estimate)
    )
epi_cor$variable <- colnames(epi_variables)

# plot the correlation to the scores
p <- epi_cor %>%
    as.data.frame() %>%
    # rownames_to_column('variable') %>%
    ggplot(aes(
        x = factor(variable, levels = epi_cor$variable), 
        y = scores,
        col = p < 0.05
        )
    ) +
    lims(y = c(0, 0.5)) +
    geom_point() +
    geom_text(aes(label = round(scores, 2)), vjust = -1) +
    theme_matt(16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # add vertical lines
    geom_vline(xintercept = 1:length(epi_variables), linetype = 'dashed', alpha = 0.2) +
    labs(x = 'Epi Variable', y = 'Correlation to \nPACE Press Score', title = 'Correlation of Epi Variables to PACE Press Score ALL SAMPLES')
p
# save the plot
ggsave(file.path(outdir, 'pace_press_scoring_epi_correlation_all.pdf'), p, width = 10, height = 10)
write.csv(epi_cor[, c('variable', 'scores', 'p')], file.path(outdir, 'pace_press_scoring_epi_correlation_all.csv'))


### DERIVATION ONLY
# map cor.test to get a p-value and correlation
subseting <- metadata$labels != 'NA'
epi_cor <- map_df(
    epi_variables[subseting, ], 
    ~cor.test(.x, scores[subseting]) %>% 
        broom::tidy() %>% 
        dplyr::select(p.value, estimate) %>% 
        dplyr::rename(p = p.value, correlation = estimate)
    )
epi_cor$variable <- colnames(epi_variables)

# plot the correlation to the scores
p <- epi_cor %>%
    as.data.frame() %>%
    # rownames_to_column('variable') %>%
    ggplot(aes(
        x = factor(variable, levels = epi_cor$variable), 
        y = correlation,
        col = p < 0.05
        )
    ) +
    lims(y = c(0, 0.5)) +
    geom_point() +
    geom_text(aes(label = round(correlation, 2)), vjust = -1) +
    theme_matt(16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # add vertical lines
    geom_vline(xintercept = 1:length(epi_variables), linetype = 'dashed', alpha = 0.2) +
    labs(x = 'Epi Variable', y = 'Correlation to \nPACE Press Score', title = 'Correlation of Epi Variables to PACE Press Score DERIVATION ONLY')
p
# save the plot
ggsave(file.path(outdir, 'pace_press_scoring_epi_correlation_derivation.pdf'), p, width = 10, height = 10)
write.csv(epi_cor[, c('variable', 'correlation', 'p')], file.path(outdir, 'pace_press_scoring_epi_correlation_derivation.csv'))



### ANTIPLATELET ONLY
# map cor.test to get a p-value and correlation
subseting <- metadata$antiplatelet_therapy == '1:Yes'
epi_cor <- map_df(
    epi_variables[subseting, ], 
    ~cor.test(.x, scores[subseting]) %>% 
        broom::tidy() %>% 
        dplyr::select(p.value, estimate) %>% 
        dplyr::rename(p = p.value, correlation = estimate)
    )
epi_cor$variable <- colnames(epi_variables)

# plot the correlation to the scores
p <- epi_cor %>%
    as.data.frame() %>%
    # rownames_to_column('variable') %>%
    ggplot(aes(
        x = factor(variable, levels = epi_cor$variable), 
        y = correlation,
        col = p < 0.05
        )
    ) +
    lims(y = c(0, 0.5)) +
    geom_point() +
    geom_text(aes(label = round(correlation, 2)), vjust = -1) +
    theme_matt(16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # add vertical lines
    geom_vline(xintercept = 1:length(epi_variables), linetype = 'dashed', alpha = 0.2) +
    labs(x = 'Epi Variable', y = 'Correlation to \nPACE Press Score', title = 'Correlation of Epi Variables to PACE Press Score ANTIPLATELET ONLY')
# save the plot
ggsave(file.path(outdir, 'pace_press_scoring_epi_correlation_on_antiplt.pdf'), p, width = 10, height = 10)
write.csv(epi_cor[, c('variable', 'correlation', 'p')], file.path(outdir, 'pace_press_scoring_epi_correlation_on_antiplt.csv'))


### ANTIPLATELET ONLY
# map cor.test to get a p-value and correlation
subseting <- metadata$antiplatelet_therapy == '2:No'
epi_cor <- map_df(
    epi_variables[subseting, ], 
    ~cor.test(.x, scores[subseting]) %>% 
        broom::tidy() %>% 
        dplyr::select(p.value, estimate) %>% 
        dplyr::rename(p = p.value, correlation = estimate)
    )
epi_cor$variable <- colnames(epi_variables)

# plot the correlation to the scores
p <- epi_cor %>%
    as.data.frame() %>%
    # rownames_to_column('variable') %>%
    ggplot(aes(
        x = factor(variable, levels = epi_cor$variable), 
        y = correlation,
        col = p < 0.05
        )
    ) +
    geom_point() +
    geom_text(aes(label = round(correlation, 2)), vjust = -1) +
    theme_matt(16) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    # add vertical lines
    geom_vline(xintercept = 1:length(epi_variables), linetype = 'dashed', alpha = 0.2) +
    labs(x = 'Epi Variable', y = 'Correlation to \nPACE Press Score', title = 'Correlation of Epi Variables to PACE Press Score OFF ANTIPLATELET ONLY')
p
# save the plot
ggsave(file.path(outdir, 'pace_press_scoring_epi_correlation_off_antiplt.pdf'), p, width = 10, height = 10)
write.csv(epi_cor[, c('variable', 'correlation', 'p')], file.path(outdir, 'pace_press_scoring_epi_correlation_off_antiplt.csv'))

#======================== HIGH/LOW PRESS SCORING LTA ========================#
# exclusions
excludes <- c('PACE062', 'PACE098', 'PACE245', 'PACE249')
stats_df <- metadata[!rownames(metadata) %in% excludes, ]

# So same breakdowns above but now as stats tables for the press scores as high / low
# add the high / low scores to the metadata
stats_df <- stats_df %>%
    mutate(
        high_low_scores = factor(ntile(scores, 2), labels = c('LOW', 'HIGH'))
        )
table(stats_df$high_low_scores)

# stats tables
# ALL SAMPLES
vars <- colnames(epi_variables)
group <- 'high_low_scores'
# make the tables
all_sample_high_low_tab <- stats_table(stats_df, group, vars, printArgs = list(nonnorm = vars))
deriv_sample_high_low_tab <- stats_table(
    stats_df %>% 
        filter(labels != 'NA') %>%
        mutate(high_low_scores = factor(ntile(scores, 2), labels = c('LOW', 'HIGH'))),
    group, vars, 
    printArgs = list(nonnorm = vars)
    )
antiplt_therapy_high_low_tab <- stats_table(
    stats_df %>% 
        filter(antiplatelet_therapy == '1:Yes') %>%
        mutate(high_low_scores = factor(ntile(scores, 2), labels = c('LOW', 'HIGH'))),
    group, vars, 
    printArgs = list(nonnorm = vars)
    )
no_antiplt_therapy_high_low_tab <- stats_table(
    stats_df %>% 
        filter(antiplatelet_therapy == '2:No') %>%
        mutate(high_low_scores = factor(ntile(scores, 2), labels = c('LOW', 'HIGH'))),
    group, vars, 
    printArgs = list(nonnorm = vars)
    )
label_high_low_tab <- stats_table(
    stats_df %>% 
        filter(labels != 'NA') %>%
        mutate(high_low_scores = factor(ntile(scores, 2), labels = c('LOW', 'HIGH'))),
    'labels', 'scores', 
    printArgs = list(nonnorm = 'scores')
    )

# save the tables
write.csv(all_sample_high_low_tab, file.path(outdir, 'pace_press_scoring_high_low_all_samples.csv'))
write.csv(deriv_sample_high_low_tab, file.path(outdir, 'pace_press_scoring_high_low_derivation.csv'))
write.csv(antiplt_therapy_high_low_tab, file.path(outdir, 'pace_press_scoring_high_low_on_antiplt.csv'))
write.csv(no_antiplt_therapy_high_low_tab, file.path(outdir, 'pace_press_scoring_high_low_off_antiplt.csv'))
write.csv(label_high_low_tab, file.path(outdir, 'pace_press_scoring_high_low_labels.csv'))

stats_df %>% dplyr::select(labels)
table(stats_df$labels)

#======================== EVENTS BY PRESS ========================#
dir.create(file.path(outdir, 'events_by_press'), showWarnings = F)
# So now let's look at the events by PRESS
censors <- grep('censor', colnames(metadata), value = T)

# make a table of the events by press
events_table_long <- metadata %>%
    dplyr::select(censors, scores) %>%
    pivot_longer(cols = c(censors), names_to = 'event', values_to = 'censor') %>%
    group_by(event, censor)
summary_table <- events_table_long %>%
    summarise(
        mean_score = mean(scores, na.rm = T),
        sd_score = sd(scores, na.rm = T),
        n = n(),
        n_censored = sum(censor == 1, na.rm = T),
        n_uncensored = sum(censor == 0, na.rm = T)
        )
write.csv(summary_table, file.path(outdir, 'events_by_press', 'summary_table.csv'))
stats_table <- events_table_long %>%
    group_by(event) %>%
    do(
        as.data.frame(stats_table(., 'censor', 'scores'))
        )
write.csv(stats_table, file.path(outdir, 'events_by_press', 'stats_table.csv'))

