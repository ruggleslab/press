###########################################################################
#
#                            duke_all_data
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-01-07
# Script Name: duke_all_data
# Output directory:
experiment <- "duke_all_data"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

# Notes about the experiment run:
notes <- "So we want to look to see how all the duke data will be in each group." #nolint

# save notes to file
write(notes, file.path(outdir, "notes.txt"))

#======================== LIBRARIES ========================#
library(lintr) #nolint
library(httpgd) #nolint
library(languageserver) #nolint
library(devtools)
library(glue)
library(magrittr)
library(ggpmisc)
library(tidyverse)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')

#======================== CODE ========================#
# Let's prep the data
duke_count <- read.csv('data/duke_validation_run3/dukerawcounttable_conv.csv', header = T, row.names = 1)
duke_meta <- read.csv('data/duke_validation_run3/dukemetatable_sel.csv', header = T, row.names = 1)
duke_long_cohort <- read.csv('data/clean/duke_longitudinal_group.csv', header = T, row.names = 1) %>% pull(.)

duke_count[1:5, 1:5]
head(duke_meta)

# load in the press genes
press <- read.csv('data/press451_genes.csv', header = T) %>% pull(.)

missing_genes <- setdiff(press, rownames(duke_count))
write.csv(missing_genes, file.path(outdir, 'duke_missing_genes.csv'))

# add the press genes to the count table
duke_count %<>% add_missing_rows(press) %>%
    filter(rownames(.) %in% press)
# make an SE now
duke_data <- make_se(duke_count, duke_meta)
duke_data <- DESeq(DESeqDataSet(duke_data, design = ~ 1))

# now let's get the normalized counts
duke_count_norm <- t(normalize_counts(duke_data, 'log2-mor'))
write.csv(duke_count_norm, file.path(outdir, 'duke_all_data.csv'))

# okay the count table is ready to go
# now we can run the run_press.py script
# this will output a file called duke_all_data.csv
# this file will have the press scores for each sample
# let's run it here
dat <- glue('{outdir}/duke_all_data.csv')
out <- glue('{outdir}/duke_all_data_scores2.csv')
system(glue('python3 code/run_press.py --data {dat} --out {out}'))

# load in the scores
duke_scores <- read.csv(out, header = TRUE, row.names = 1)
head(duke_scores)
skimr::skim(duke_scores)
#======================== Break down Scores by Metadata ========================#
vars <- c('characteristic__subject_id', 'characteristic__treatment_drug_label', 'characteristic__visit_number' , 'long_cohort')

# add the scores to the metadata
meta <- cbind(duke_meta, duke_scores)
meta %<>% mutate(
    hyper_responder = factor(
        case_when(
            characteristic__epi_max_05 > 60 ~ 'hyper',
            characteristic__epi_max_05 < 40 ~ 'normal',
            TRUE ~ NA_character_
        ), levels = c('normal', 'hyper')),
    long_cohort = factor(case_when(
            rownames(.) %in% duke_long_cohort ~ hyper_responder,
            TRUE ~ NA_character_
            ), levels = c('normal', 'hyper')),
    characteristic__treatment_drug_label = factor(case_when(
            characteristic__visit_number == 1 ~ 'Baseline',
            characteristic__visit_number == 4 ~ 'Washout',
            TRUE ~ characteristic__treatment_drug_label),
        levels = c('Baseline', '81 mg Aspirin', '325 mg Aspirin', '180 mg Ticagrelor', 'Washout'))
    )
# meta %>%
#     filter(characteristic__subject_id %in% ids & hyper_baseline == 'hyper') %>%
#     rstatix::t_test(scores ~ characteristic__treatment_drug_label, paired = TRUE, detailed = TRUE)

buckets <- meta %>% 
    filter(characteristic__visit_number == 1) %>% 
    dplyr::select(scores, characteristic__subject_id) %>%
    mutate(press_bucket = ntile(scores, 4))
# histogram of the buckets
histogram <- ggplot(buckets, aes(x = scores, fill = factor(press_bucket))) +
    geom_histogram() +
    theme_bw(20) +
    labs(x = 'PRESS Score', y = 'Count', title = 'Histogram of PRESS Scores') +
    theme(legend.position = 'none')
ggsave(file.path(outdir, 'baseline_histogram.png'), histogram, width = 8, height = 6)
# match the buckets to the correct subject
meta %<>% left_join(dplyr::select(buckets, -scores), by = 'characteristic__subject_id')
colnames(meta)

# make a histogram of the scores and color by the bucket
histogram <- ggplot(meta, aes(x = scores, fill = factor(press_bucket))) +
    geom_histogram() +
    theme_bw(20) +
    labs(x = 'PRESS Score', y = 'Count', title = 'Histogram of PRESS Scores') +
    theme(legend.position = 'none')
ggsave(file.path(outdir, 'histogram.png'), histogram, width = 8, height = 6)

ids <- meta %>% group_by(characteristic__subject_id) %>% summarise(n = n()) %>% filter(n == 5) %>% pull(characteristic__subject_id)


# get the hyper responders for baseline
hyper_baseline <- meta %>% 
    filter(characteristic__visit_number == 1) %>%
    dplyr::select(characteristic__subject_id, hyper_responder)
# make a column of the baseline hyper responders by matching the subject id
meta %<>% mutate(
    hyper_baseline = factor(hyper_baseline[match(characteristic__subject_id, hyper_baseline$characteristic__subject_id), 'hyper_responder'],
                            levels = c('normal', 'hyper'))
    )
# meta %<>% mutate_at(all_of(vars), as.character)


drug_label_plot <- ggplot(meta, aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    theme(legend.position = 'none') +
    stat_compare_means(
        method = 't.test',
        comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
        size = 6
        )
ggsave(file.path(outdir, 'drug_label_plot.png'), drug_label_plot, width = 10, height = 5)

drug_label_plot_paired <- ggplot(
    meta %>% filter(characteristic__subject_id %in% ids), 
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id, color = factor(press_bucket)), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(
        method = 't.test',
        comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
        size = 6,
        paired = TRUE
        )
ggsave(file.path(outdir, 'drug_label_plot_paired.png'), drug_label_plot_paired, width = 10, height = 5)

drug_label_plot_paired_hyper <- ggplot(
    meta %>% filter(characteristic__subject_id %in% ids & hyper_baseline == 'hyper'), 
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(
        method = 't.test',
        comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
        size = 6,
        paired = TRUE
        )
ggsave(file.path(outdir, 'drug_label_plot_paired_hyper.png'), drug_label_plot_paired_hyper, width = 10, height = 5)

drug_label_plot_paired_bucket <- ggplot(
    meta %>% filter(press_bucket == '4' & characteristic__subject_id %in% ids),
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(
        method = 't.test',
        comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
        size = 6,
        paired = TRUE
        )
ggsave(file.path(outdir, 'drug_label_plot_paired_bucket.png'), drug_label_plot_paired_bucket, width = 10, height = 5)

visit_num_plot <- ggplot(meta, aes(x = characteristic__visit_number, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    theme_bw(16) +
    labs(x = 'Visit Number', y = 'PRESS Score', title = 'PRESS Scores by Visit Number') +
    theme(legend.position = 'right')
ggsave(file.path(outdir, 'visit_num_plot.png'), visit_num_plot, width = 10, height = 5)

ids <- meta %>% group_by(characteristic__subject_id) %>% summarise(n = n()) %>% filter(n == 5) %>% pull(characteristic__subject_id)
drug_hyper_plot <- meta %>%
    drop_na(hyper_baseline) %>%
        ggplot(aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
        geom_boxplot() +
        geom_jitter(width = 0.1, alpha = 0.5, aes(color = hyper_baseline)) +
        theme_bw(20) +
        labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
        theme(legend.position = 'bottom') +
        # facet_grid(hyper_baseline ~ .) +
        stat_compare_means(
            method = 't.test',
            paired = TRUE,
            comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
            size = 6
        )
ggsave(file.path(outdir, 'drug_hyper_plot.png'), drug_hyper_plot, width = 12, height = 10)
drug_hyper_tab <- meta %>%
    drop_na(hyper_baseline) %>%
    group_by(characteristic__treatment_drug_label, hyper_baseline) %>%
    summarise(n = n(), .groups = 'drop')
write.csv(drug_hyper_tab, file.path(outdir, 'drug_hyper_tab.csv'))

drug_hyper_hyper_plot <- meta %>% 
    filter(hyper_baseline == 'hyper') %>%
    ggplot(aes(x = characteristic__treatment_drug_label, y = scores, fill = hyper_baseline)) +
    geom_boxplot() +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    theme(legend.position = 'top') +
    stat_compare_means(
        method = 't.test',
        comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
        size = 6
    )
ggsave(file.path(outdir, 'drug_hyper_hyper_plot.png'), drug_hyper_hyper_plot, width = 8, height = 5)

drug_hyper_normal_plot <- meta %>% 
    filter(hyper_baseline == 'normal') %>%
    ggplot(aes(x = characteristic__treatment_drug_label, y = scores, fill = hyper_baseline)) +
    geom_boxplot() +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    theme(legend.position = 'top') +
    stat_compare_means(
        method = 't.test',
        comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)),
        size = 6
    )
ggsave(file.path(outdir, 'drug_hyper_normal_plot.png'), drug_hyper_normal_plot, width = 8, height = 5)

platelet_function_plot <- meta %>%
    ggplot(aes(x = characteristic__PlateletFunctionScore, y = scores, color = characteristic__treatment_drug_label)) +
        geom_point() +
        theme_bw(24) +
        theme(legend.position = 'bottom') +
        labs(x = 'Platelet Function Score', y = 'PRESS Score', title = 'PRESS Scores by Platelet Function Score') +
        stat_smooth(method = 'lm', se = TRUE) +
        stat_cor(method = 'pearson', size = 6, show.legend = FALSE)
ggsave(file.path(outdir, 'platelet_function_plot.png'), platelet_function_plot, width = 8, height = 6)

all_hyper_hypo_plot <- meta %>% 
    drop_na(hyper_baseline) %>%
    filter(characteristic__treatment_drug_label == 'none') %>%
    ggplot(aes(x = hyper_baseline, y = scores, fill = hyper_baseline)) +
    geom_boxplot() +
    theme_bw(24) +
    theme(legend.position = 'none') +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    stat_compare_means(method = 't.test', size = 6)
ggsave(file.path(outdir, 'all_hyper_hypo_plot.png'), all_hyper_hypo_plot, width = 6, height = 6)

all_hyper_hypo_plot <- meta %>% 
    drop_na(cohort) %>%
    ggplot(aes(x = hyper_baseline, y = scores, fill = cohort)) +
    geom_boxplot() +
    theme_bw(24) +
    theme(legend.position = 'bottom') +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    stat_compare_means(method = 't.test', size = 6)
ggsave(file.path(outdir, 'cohort_drug_plot.png'), all_hyper_hypo_plot, width = 6, height = 6)

long_hyper_hypo_plot <- meta %>% 
    drop_na(long_cohort) %>%
    filter(characteristic__treatment_drug_label == 'none') %>%
    ggplot(aes(x = long_cohort, y = scores, fill = long_cohort)) +
    geom_boxplot() +
    theme_bw(24) +
    theme(legend.position = 'none') +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    stat_compare_means(method = 't.test', size = 6)
ggsave(file.path(outdir, 'long_hyper_hypo_plot.png'), long_hyper_hypo_plot, width = 6, height = 6)

#======================== ROC ========================#
# get the auc to double check if it's all good
auc <- with(meta, pROC::roc(long_cohort, preds))
auc
with(meta, table(hyper_baseline, cohort))

# make a stats table for duke data
duke_stats <- stats_table(
    meta %>% dplyr::select(-c(characteristic__subject_id, characteristic__treatment_drug)), 
    'cohort'
    )
write.csv(duke_stats, file.path(outdir, 'duke_stats.csv'))

duke_stats_by_drug <- stats_table(
    meta %>% dplyr::select(-c(characteristic__subject_id, characteristic__treatment_drug, long_cohort)), 
    'characteristic__treatment_drug_label'
    )
write.csv(duke_stats_by_drug, file.path(outdir, 'duke_stats_by_drug.csv'))

#======================== PRESS Stability in DUKE ========================#
# Okay another cool question to ask -- how stable is press at baseline versus washout
# this is group1 versus group2 and might as well fill by baseline hyper responder
# pivot data wider to have a group1 and group2 column of the scores
ids <- meta[duke_long_cohort, 'characteristic__subject_id']
press_stability <- meta %>% filter(cohort %in% c('group1', 'group2')) %>%
    drop_na(hyper_baseline) %>%
    # filter(characteristic__subject_id %in% ids) %>%
    pivot_wider(id_cols = c(hyper_baseline, characteristic__subject_id), names_from = cohort, values_from = scores) %>%
    dplyr::select(group1, group2, hyper_baseline, characteristic__subject_id)
with(press_stability, table(hyper_baseline))

press_g1_g2_plot <- ggplot(press_stability, aes(x = group1, y = group2)) +
    geom_point() +
    theme_bw(24) +
    theme(legend.position = 'bottom') +
    labs(x = 'Group 1', y = 'Group 2', title = 'PRESS Scores by Group') +
    stat_smooth(method = 'lm', se = TRUE) +
    stat_poly_eq(use_label(c("R2", "R2.CI", "P", "method")), formula = y ~ x, parse = TRUE, vjust = 3, size = 6) +
    stat_cor(method = 'spearman', size = 6, show.legend = FALSE)
press_g1_g2_plot
ggsave(file.path(outdir, 'press_g1_g2_plot.png'), press_g1_g2_plot, width = 8, height = 6)


#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
