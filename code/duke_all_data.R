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

#======================== LIBRARIES ========================
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

#======================== CODE ========================
# Let's prep the data
duke_count <- read.csv('data/duke_validation_run3/dukerawcounttable_conv.csv', header = T, row.names = 1)
duke_meta <- read.csv('data/GSE158765/MetaData_samples_withreactivitytraits.csv', header = T, row.names = 1)
duke_long_cohort <- read.csv('data/clean/duke_longitudinal_group.csv', header = T, row.names = 1) %>% pull(.)

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
#======================== Break down Scores by Metadata ========================
dir.create(file.path(outdir, 'breakdowns_metadata'), showWarnings = F)
vars <- c('characteristic__subject_id', 'characteristic__treatment_drug_label', 'characteristic__visit_number' , 'long_cohort')

# add the scores to the metadata
meta <- cbind(as.data.frame(colData(duke_data)), duke_scores)
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
        levels = c('Baseline', '81 mg Aspirin', '325 mg Aspirin', '180 mg Ticagrelor', 'Washout')),
    cohort = case_when(
        characteristic__visit_number == 1 ~ 'group1',
        characteristic__visit_number == 4 ~ 'group2',
        characteristic__treatment_drug == 1 ~ 'group3',
        characteristic__treatment_drug == 2 ~ 'group4',
        characteristic__treatment_drug == 3 ~ 'group5'
    )
    )
# meta %>%
#     filter(characteristic__subject_id %in% ids & hyper_baseline == 'hyper') %>%
#     rstatix::t_test(scores ~ characteristic__treatment_drug_label, paired = TRUE, detailed = TRUE)

buckets <- meta %>% 
    filter(characteristic__visit_number == 1) %>% 
    dplyr::select(scores, characteristic__subject_id) %>%
    mutate(press_bucket = ntile(scores, 2))
# histogram of the buckets
histogram <- ggplot(buckets, aes(x = scores, fill = factor(press_bucket))) +
    geom_histogram() +
    theme_bw(20) +
    labs(x = 'PRESS Score', y = 'Count', title = 'Histogram of PRESS Scores') +
    theme(legend.position = 'none')
ggsave(file.path(outdir, 'breakdowns_metadata', 'baseline_histogram.png'), histogram, width = 8, height = 6)
# match the buckets to the correct subject
rows <- rownames(meta)
meta %<>% left_join(dplyr::select(buckets, -scores), by = 'characteristic__subject_id')
rownames(meta) <- rows
colnames(meta)

# make a histogram of the scores and color by the bucket
histogram <- ggplot(meta, aes(x = scores, fill = factor(press_bucket))) +
    geom_histogram() +
    theme_bw(20) +
    labs(x = 'PRESS Score', y = 'Count', title = 'Histogram of PRESS Scores') +
    theme(legend.position = 'none')
ggsave(file.path(outdir, 'breakdowns_metadata', 'histogram.png'), histogram, width = 8, height = 6)

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
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot.png'), drug_label_plot, width = 10, height = 5)

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
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired.png'), drug_label_plot_paired, width = 10, height = 5)

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
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired_hyper.png'), drug_label_plot_paired_hyper, width = 10, height = 5)

drug_label_plot_paired_bucket <- ggplot(
    meta %>% filter(press_bucket == '2' & characteristic__subject_id %in% ids),
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(method = 't.test', size = 6, paired = TRUE, comparisons = pairwise_combos(as.character(meta$characteristic__treatment_drug_label)))
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired_bucket.png'), drug_label_plot_paired_bucket, width = 10, height = 5)

drug_label_plot_paired_bucket_81asa <- ggplot(
    meta %>% filter(press_bucket == '2' & characteristic__subject_id %in% ids & characteristic__treatment_drug_label %in% c('81 mg Aspirin', 'Baseline')),
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(method = 't.test', size = 6, paired = TRUE)
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired_bucket_81asa.png'), drug_label_plot_paired_bucket_81asa, width = 5, height = 5)

drug_label_plot_paired_bucket_325asa <- ggplot(
    meta %>% filter(press_bucket == '2' & characteristic__subject_id %in% ids & characteristic__treatment_drug_label %in% c('325 mg Aspirin', 'Baseline')),
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(method = 't.test', size = 6, paired = TRUE)
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired_bucket_325asa.png'), drug_label_plot_paired_bucket_325asa, width = 5, height = 5)

drug_label_plot_paired_bucket_180tic <- ggplot(
    meta %>% filter(press_bucket == '2' & characteristic__subject_id %in% ids & characteristic__treatment_drug_label %in% c('180 mg Ticagrelor', 'Baseline')),
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(method = 't.test', size = 6, paired = TRUE)
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired_bucket_180tic.png'), drug_label_plot_paired_bucket_180tic, width = 5, height = 5)
with(meta %>% filter(cohort == 'group1'), table(press_bucket, hyper_baseline))

drug_label_plot_paired_bucket_washout <- ggplot(
    meta %>% filter(press_bucket == '2' & characteristic__subject_id %in% ids & characteristic__treatment_drug_label %in% c('Washout', '180 mg Ticagrelor')),
    aes(x = characteristic__treatment_drug_label, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
    theme_bw(20) +
    labs(x = 'Group', y = 'PRESS Score', title = NULL) +
    theme(legend.position = 'none') +
    stat_compare_means(method = 't.test', size = 6, paired = TRUE)
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_label_plot_paired_bucket_washout.png'), drug_label_plot_paired_bucket_washout, width = 5, height = 5)

visit_num_plot <- ggplot(meta, aes(x = characteristic__visit_number, y = scores, fill = characteristic__treatment_drug_label)) +
    geom_boxplot() +
    theme_bw(16) +
    labs(x = 'Visit Number', y = 'PRESS Score', title = 'PRESS Scores by Visit Number') +
    theme(legend.position = 'right')
ggsave(file.path(outdir, 'breakdowns_metadata', 'visit_num_plot.png'), visit_num_plot, width = 10, height = 5)

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
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_hyper_plot.png'), drug_hyper_plot, width = 12, height = 10)
drug_hyper_tab <- meta %>%
    drop_na(hyper_baseline) %>%
    group_by(characteristic__treatment_drug_label, hyper_baseline) %>%
    summarise(n = n(), .groups = 'drop')
write.csv(drug_hyper_tab, file.path(outdir, 'breakdowns_metadata', 'drug_hyper_tab.csv'))

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
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_hyper_hyper_plot.png'), drug_hyper_hyper_plot, width = 8, height = 5)

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
ggsave(file.path(outdir, 'breakdowns_metadata', 'drug_hyper_normal_plot.png'), drug_hyper_normal_plot, width = 8, height = 5)

platelet_function_plot <- meta %>%
    ggplot(aes(x = characteristic__PlateletFunctionScore, y = scores, color = characteristic__treatment_drug_label)) +
        geom_point() +
        theme_bw(24) +
        theme(legend.position = 'bottom') +
        labs(x = 'Platelet Function Score', y = 'PRESS Score', title = 'PRESS Scores by Platelet Function Score') +
        stat_smooth(method = 'lm', se = TRUE) +
        stat_cor(method = 'pearson', size = 6, show.legend = FALSE)
ggsave(file.path(outdir, 'breakdowns_metadata', 'platelet_function_plot.png'), platelet_function_plot, width = 8, height = 6)

all_hyper_hypo_plot <- meta %>% 
    drop_na(hyper_baseline) %>%
    filter(characteristic__treatment_drug_label == 'none') %>%
    ggplot(aes(x = hyper_baseline, y = scores, fill = hyper_baseline)) +
    geom_boxplot() +
    theme_bw(24) +
    theme(legend.position = 'none') +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    stat_compare_means(method = 't.test', size = 6)
ggsave(file.path(outdir, 'breakdowns_metadata', 'all_hyper_hypo_plot.png'), all_hyper_hypo_plot, width = 6, height = 6)

all_hyper_hypo_plot <- meta %>% 
    drop_na(cohort) %>%
    ggplot(aes(x = hyper_baseline, y = scores, fill = cohort)) +
    geom_boxplot() +
    theme_bw(24) +
    theme(legend.position = 'bottom') +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    stat_compare_means(method = 't.test', size = 6)
ggsave(file.path(outdir, 'breakdowns_metadata', 'cohort_drug_plot.png'), all_hyper_hypo_plot, width = 6, height = 6)

long_hyper_hypo_plot <- meta %>% 
    drop_na(long_cohort) %>%
    filter(characteristic__treatment_drug_label == 'none') %>%
    ggplot(aes(x = long_cohort, y = scores, fill = long_cohort)) +
    geom_boxplot() +
    theme_bw(24) +
    theme(legend.position = 'none') +
    labs(x = 'Group', y = 'PRESS Score', title = 'PRESS Scores by Group') +
    stat_compare_means(method = 't.test', size = 6)
ggsave(file.path(outdir, 'breakdowns_metadata', 'long_hyper_hypo_plot.png'), long_hyper_hypo_plot, width = 6, height = 6)

#======================== Tx ON PRESS ========================
dir.create(file.path(outdir, 'tx'), showWarnings = F)
# Effect of ASA on PRESS in DUKE (low & high dose)
# - overall
# - Pts with PRESS > median
# - Pts with hyper prediction
# Effect on:
# - PRESS Score
# - Platelet Score
# - AA LTA
# - Epi LTA
# - ADP LTA
# - PLT count
# - MPV

tx_groups <- c('81 mg Aspirin', '325 mg Aspirin', '180 mg Ticagrelor')
control_group <- c('Baseline', 'Washout')
pop_groups <- c('all', 'hyper', 'normal', 'above_median', 'below_median')
all_samples <- meta %>% pull(characteristic__subject_id)
hyper_samples <- meta %>% 
    filter(cohort == 'group1') %>%
    filter(scores > 0.38) %>% 
    pull(characteristic__subject_id)
normal_samples <- meta %>% 
    filter(cohort == 'group1') %>%
    filter(scores < 0.38) %>% 
    pull(characteristic__subject_id)
above_median_samples <- meta %>% filter(press_bucket == 2) %>% pull(characteristic__subject_id)
below_median_samples <- meta %>% filter(press_bucket == 1) %>% pull(characteristic__subject_id)
# let's perform the mediation analysis on a matched set of features
regex <- 'scores|PlateletFunctionScore|__(aa|arach|epi|adp|coll)|abc_(plt|mpv)'
plt_vars <- grep(regex, colnames(meta), value = T, perl = T)
plt_vars <- plt_vars[!grepl("trt$", plt_vars)]

for (pop in pop_groups) {
    dir.create(file.path(outdir, 'tx', pop), showWarnings = F, recursive = T)
    overall_res <- data.frame()
    for (tx in tx_groups) {
        for (ctrl in control_group) {
            df <- meta %>% 
                filter(characteristic__treatment_drug_label %in% c(tx, ctrl)) %>%
                arrange(characteristic__subject_id, characteristic__treatment_drug_label)
            paired_samples <- df %>% group_by(characteristic__subject_id) %>% filter(n() == 2) %>% pull(characteristic__subject_id)
            df %<>% 
                filter(characteristic__subject_id %in% paired_samples) %>%
                filter(characteristic__subject_id %in% switch(pop,
                    'all' = all_samples,
                    'hyper' = hyper_samples,
                    'normal' = normal_samples,
                    'above_median' = above_median_samples,
                    'below_median' = below_median_samples
                ))
            write.csv(df, file.path(outdir, 'tx', pop, glue('{tx}_{ctrl}.csv')))
            dir.create(file.path(outdir, 'tx', pop, tx), showWarnings = F)
            plt_vars_df <- map_df(plt_vars, ~{
                var <- .x
                paired_samples <- df %>% 
                    drop_na(!!sym(var)) %>%
                    group_by(characteristic__subject_id) %>% filter(n() == 2) %>% 
                    pull(characteristic__subject_id)

                paired_label_plot <- ggplot(
                    df %>% filter(characteristic__subject_id %in% paired_samples), 
                    aes(
                        x = characteristic__treatment_drug_label, 
                        y = !!sym(var),
                        color = characteristic__treatment_drug_label
                        )
                    ) +
                    geom_boxplot() +
                    geom_jitter(width = 0.01) +
                    geom_line(aes(group = characteristic__subject_id), alpha = 0.5) +
                    theme_bw(20) +
                    labs(x = NULL, y = toupper(gsub('characteristic__', '', var)), title = NULL) +
                    theme(legend.position = 'none') +
                    stat_compare_means(
                        method = 't.test',
                        size = 6,
                        paired = TRUE
                    )
                ggsave(file.path(outdir, 'tx', pop, tx, glue('{tx}_{ctrl}_{var}.png')), paired_label_plot, width = 8, height = 6)

                t_test_res <- df %>% 
                    filter(characteristic__subject_id %in% paired_samples) %>%
                    mutate(characteristic__treatment_drug_label = as.character(characteristic__treatment_drug_label)) %>%
                    rstatix::t_test(as.formula(glue('{var} ~ characteristic__treatment_drug_label')), paired = TRUE, detailed = TRUE) %>%
                    mutate(var = gsub('characteristic__', '', var)) %>%
                    dplyr::select(
                        var, group1, group2, n = n1, 
                        estimate, p, conf.low, conf.high, 
                        method, alternative
                        )
                return(t_test_res)
            })
            overall_res <- rbind(overall_res, plt_vars_df)
            write.csv(plt_vars_df, file.path(outdir, 'tx', pop, glue('{tx}_{ctrl}_stats.csv')))
            sumPlot <- ggplot(plt_vars_df, aes(x = var, y = estimate, color = p < 0.05)) +
                geom_point() +
                geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
                geom_hline(yintercept = 0, linetype = 'dashed') +
                coord_flip() +
                theme_bw(20) +
                labs(x = NULL, y = 'Mean Difference', title = 'Effect of Treatment on PRESS') +
                theme(legend.position = 'top')
            ggsave(file.path(outdir, 'tx', pop, glue('{tx}_{ctrl}_sumplot.png')), sumPlot, width = 8, height = 8)
        }
    }
    write.csv(overall_res, file.path(outdir, 'tx', pop, 'overall_stats.csv'))
    tmp <- overall_res %>% 
        dplyr::select(var, group1, group2, p) %>%
        mutate(p = -log10(p)) %>%
        pivot_wider(names_from = c(group1, group2), values_from = p) %>%
        arrange(var) %>%
        column_to_rownames('var') %>%
        drop_na() %>%
        t()
    hm <- Heatmap(
        tmp,
        name = '-log10(p)',
        col = circlize::colorRamp2(c(0, 5), c('white', 'red')),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "black", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(ifelse(tmp[i, j] > -log10(0.05), '*', ' '), x, y, gp = gpar(fontsize = 20))
        }
    )
    pdf(file.path(outdir, 'tx', pop, 'overall_heatmap.pdf'), width = 12, height = 6)
    draw(hm)
    dev.off()
}


#======================== ROC ========================
# # get the auc to double check if it's all good
# auc <- with(meta, pROC::roc(long_cohort, preds))
# auc
# with(meta, table(hyper_baseline, cohort))

# # make a stats table for duke data
# duke_stats <- stats_table(
#     meta %>% dplyr::select(-c(characteristic__subject_id, characteristic__treatment_drug)), 
#     'cohort'
#     )
# write.csv(duke_stats, file.path(outdir, 'duke_stats.csv'))

# duke_stats_by_drug <- stats_table(
#     meta %>% dplyr::select(-c(characteristic__subject_id, characteristic__treatment_drug, long_cohort)), 
#     'characteristic__treatment_drug_label'
#     )
# write.csv(duke_stats_by_drug, file.path(outdir, 'duke_stats_by_drug.csv'))

#======================== PRESS Stability in DUKE ========================
dir.create(file.path(outdir, 'stability'), showWarnings = F)
# Okay another cool question to ask -- how stable is press at baseline versus washout
# this is group1 versus group2 and might as well fill by baseline hyper responder
# pivot data wider to have a group1 and group2 column of the scores
ids <- meta[duke_long_cohort, 'characteristic__subject_id']
press_stability <- meta %>% filter(cohort %in% c('group1', 'group2')) %>%
    drop_na(hyper_baseline) %>%
    # filter(characteristic__subject_id %in% ids) %>%
    pivot_wider(id_cols = c(hyper_baseline, characteristic__subject_id), names_from = cohort, values_from = scores) %>%
    dplyr::select(group1, group2, hyper_baseline, characteristic__subject_id) %>%
    mutate(concordant = case_when(
        group1 > 0.38 & group2 > 0.38 ~ 'concordant',
        group1 < 0.38 & group2 < 0.38 ~ 'concordant',
        TRUE ~ 'discordant'
    ))
write.csv(press_stability, file.path(outdir, 'stability', 'press_stability.csv'))

press_g1_g2_plot <- ggplot(press_stability, aes(x = group1, y = group2)) +
    geom_point(aes(color = concordant)) +
    theme_bw(24) +
    theme(legend.position = 'bottom') +
    geom_vline(xintercept = 0.38, linetype = 'dashed') +
    geom_hline(yintercept = 0.38, linetype = 'dashed') +
    labs(x = 'Group 1', y = 'Group 2', title = 'PRESS Scores by Group') +
    stat_smooth(method = 'lm', se = TRUE) +
    stat_poly_eq(use_label(c("R2", "R2.CI", "P", "method")), formula = y ~ x, parse = TRUE, vjust = 3, size = 6) +
    stat_cor(method = 'spearman', size = 6, show.legend = FALSE)
ggsave(file.path(outdir, 'stability', 'press_g1_g2_plot.png'), press_g1_g2_plot, width = 8, height = 6)

#======================== PRESS Change from Therapy ========================
dir.create(file.path(outdir, 'dTx'), showWarnings = F)
dTx_data <- meta %>% 
    pivot_wider(
        id_cols = c(hyper_baseline, characteristic__subject_id, press_bucket), 
        names_from = characteristic__treatment_drug_label, 
        values_from = scores
        ) %>%
    mutate(
        `81 mg Aspirin - Baseline` = `81 mg Aspirin` - Baseline,
        `325 mg Aspirin - Baseline` = `325 mg Aspirin` - Baseline,
        `180 mg Ticagrelor - Baseline` = `180 mg Ticagrelor` - Baseline,
        `Washout - Baseline` = Washout - Baseline
    ) %>%
    pivot_longer(
        cols = c(`81 mg Aspirin - Baseline`, `325 mg Aspirin - Baseline`, `180 mg Ticagrelor - Baseline`, `Washout - Baseline`),
        names_to = 'therapy',
        values_to = 'change'
    )
write.csv(dTx_data, file.path(outdir, 'dTx', 'dTx_data.csv'))

dTx_plot <- ggplot(dTx_data, aes(x = change, y = therapy)) +
    geom_boxplot() +
    theme_bw(24) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    theme(legend.position = 'bottom') +
    labs(x = 'Change in PRESS Score', y = NULL, title = 'PRESS Score Change by Group')
ggsave(file.path(outdir, 'dTx', 'dTx_plot.png'), dTx_plot, width = 8, height = 6)

dTx_hyper_plot <- dTx_data %>% 
    ggplot(aes(x = change, y = therapy, color = hyper_baseline)) +
    geom_boxplot() +
    theme_bw(24) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    theme(legend.position = 'bottom') +
    labs(x = 'Change in PRESS Score', y = NULL, title = 'PRESS Score Change by Group')
ggsave(file.path(outdir, 'dTx', 'dTx_hyper_plot.png'), dTx_hyper_plot, width = 8, height = 6)

dTx_bucket_plot <- dTx_data %>% 
    ggplot(aes(x = change, y = therapy, color = factor(press_bucket))) +
    geom_boxplot() +
    theme_bw(24) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    theme(legend.position = 'bottom') +
    labs(x = 'Change in PRESS Score', y = NULL, title = 'PRESS Score Change by Group', color = 'PRESS Bucket')
ggsave(file.path(outdir, 'dTx', 'dTx_bucket_plot.png'), dTx_bucket_plot, width = 8, height = 6)


#======================== Platelet Variables ========================
dir.create(file.path(outdir, 'platelet_variables'), showWarnings = F)
# let's compare press to the other platelet variables that are included in the duke data
colnames(meta)

regex <- 'PlateletFunctionScore|__(aa|arach|epi|adp|coll)|abc_(plt|mpv)'
variables <- grep(regex, colnames(meta), value = TRUE)
for (var in variables) {
    baseline_plot <- ggplot(
        meta %>% drop_na(hyper_baseline) %>% filter(cohort == 'group1'), 
        aes(x = .data[[var]], y = scores)
        ) +
        geom_point() +
        theme_bw(24) +
        theme(legend.position = 'none') +
        labs(x = var, y = 'PRESS Score', title = 'PRESS Scores by Platelet Function Score') +
        stat_smooth(method = 'lm', se = TRUE) +
        stat_cor(method = 'spearman', size = 6)
    ggsave(file.path(outdir, 'platelet_variables', glue('{var}_baseline_plot.png')), baseline_plot, width = 8, height = 6)

    follow_up_plot <- ggplot(
        meta %>% drop_na(hyper_baseline) %>% filter(cohort == 'group2'), 
        aes(x = .data[[var]], y = scores)
        ) +
        geom_point() +
        theme_bw(24) +
        theme(legend.position = 'none') +
        labs(x = var, y = 'PRESS Score', title = 'PRESS Scores by Platelet Function Score') +
        stat_smooth(method = 'lm', se = TRUE) +
        stat_cor(method = 'spearman', size = 6)
    ggsave(file.path(outdir, 'platelet_variables', glue('{var}_follow_up_plot.png')), follow_up_plot, width = 8, height = 6)
}
or_table_baseline <- purrr::map_df(
    variables[-1], function(y) {
        f <- as.formula(glue('{y} ~ scores'))
        dat <- meta %>% drop_na(hyper_baseline) %>% filter(cohort == 'group1') %>% mutate_at(vars(y), scale)
        tryCatch({
            m <- glm(f, data = dat)
            broom::tidy(m) %>% 
                filter(term == 'scores') %>%
                    mutate(term = y, conf.up = estimate + 1.96 * std.error, conf.low = estimate - 1.96 * std.error) %>%
                    dplyr::select(term, estimate, conf.low, conf.up, p = p.value)
        }, error = function(e) {
            return(data.frame(term = y, estimate = 0, std.error = 0, p.value = 1))
        })
    })
write.csv(or_table_baseline, file.path(outdir, 'platelet_variables', 'or_table_baseline.csv'))
p <- ggplot(or_table_baseline %>% filter(!grepl('trt', term)), aes(x = estimate, y = term, color = p < 0.05)) +
    geom_point() +
    geom_linerange(aes(xmin = conf.low, xmax = conf.up)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    lims(x = c(-0.5, 0.5)) +
    theme_bw(12) +
    theme(legend.position = 'bottom')
ggsave(file.path(outdir, 'platelet_variables', 'or_table_baseline_plot.png'), p, width = 8, height = 6)
or_table_follow_up <- purrr::map_df(
    variables, function(y) {
        f <- as.formula(glue('{y} ~ scores'))
        dat <- meta %>% drop_na(hyper_baseline) %>% filter(cohort == 'group2') %>% mutate_at(vars(y), scale)
        tryCatch({
            m <- glm(f, data = dat)
            broom::tidy(m) %>% 
                filter(term == 'scores') %>%
                    mutate(term = y, conf.up = estimate + 1.96 * std.error, conf.low = estimate - 1.96 * std.error) %>%
                    dplyr::select(term, estimate, conf.low, conf.up, p = p.value)
        }, error = function(e) {
            return(data.frame(term = y, estimate = 0, std.error = 0, p.value = 1))
        })
    })
write.csv(or_table_follow_up, file.path(outdir, 'platelet_variables', 'or_table_follow_up.csv'))
p <- ggplot(or_table_follow_up %>% filter(!grepl('trt', term)), aes(x = estimate, y = term, color = p < 0.05)) +
    geom_point() +
    geom_linerange(aes(xmin = conf.low, xmax = conf.up)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    lims(x = c(-0.5, 0.5)) +
    theme_bw(12) +
    theme(legend.position = 'bottom')
ggsave(file.path(outdir, 'platelet_variables', 'or_table_follow_up_plot.png'), p, width = 8, height = 6)

#======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))




#======================== Scratch Code ========================
ggplot(meta %>% drop_na(long_cohort), aes(x = long_cohort, y = scores)) +
    geom_boxplot() +
    stat_compare_means(method = 't.test') +
    lims(y = c(-2, 3))
