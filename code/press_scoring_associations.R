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

#======================== LIBRARIES ========================
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

#======================== PACE ========================
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

# add on the bleeding metadata
metadata_addon <- read.csv('data/metadata_9_20220422.csv', row.names = 1)
rownames(metadata_addon) <- gsub('_', '', rownames(metadata_addon))
columns_to_add <- grep('bleed', colnames(metadata_addon), value = T)
metadata[rownames(metadata_addon), columns_to_add] <- metadata_addon[ , columns_to_add]

# get the intersection of the two
matched_samples <- intersect(rownames(press_scores), rownames(metadata))

# add the scores to the metadata
metadata_all <- cbind(metadata[matched_samples,], press_scores[matched_samples,])

# add 3 and 4 tiles based on the scores to the metadata
metadata <- metadata_all %>% 
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

#======================== PRESS SCORING LTA & FLOW ========================
dir.create(file.path(outdir, 'plt_measures'), showWarnings = F)
# exclusions
excludes <- c('PACE062', 'PACE098', 'PACE245', 'PACE249')
stats_df <- metadata[!rownames(metadata) %in% excludes, ]
stats_df <- stats_df %>% 
    mutate(
        high_low_scores = factor(ntile(preds, 2), labels = c('LOW', 'HIGH')),
        labels = ifelse(labels == "NA", NA, labels)
        )
write.csv(stats_df, file.path(outdir, 'plt_measures', 'pace_metadata_with_press_scoring.csv'))

group <- 'high_low_scores'
lta <- grep('^(epi_|col_|adp_|ser_|aa_*)(?!.*slope.*|.*lag.*|.*mfi.*|.*mlty.*)', colnames(metadata), value = T, perl = T)
flow <- grep('(mfi_n|gate_n)$', colnames(metadata), value = T, perl = T)
varsList <- list(lta = lta, flow = flow)
for (vname in names(varsList)) {
    dir.create(file.path(outdir, 'plt_measures', vname), showWarnings = F)
    vars <- varsList[[vname]]
    all_sample_high_low_tab <- stats_table(stats_df, group, vars, printArgs = list(nonnorm = vars))
    write.csv(all_sample_high_low_tab, file.path(outdir, 'plt_measures', vname, 'pace_press_scoring_high_low_all_samples.csv'))
    deriv_sample_high_low_tab <- stats_table(filter(stats_df, labels != 'NA'), group, vars, printArgs = list(nonnorm = vars))
    write.csv(deriv_sample_high_low_tab, file.path(outdir, 'plt_measures', vname, 'pace_press_scoring_high_low_deriv_samples.csv'))
    # Now let's look at it continuously with all samples
    odds_ratio_table <- purrr::map_df(
        vars, 
        function(x) {
            f <- as.formula(paste0(x, ' ~ scores'))
            dat <- stats_df %>% drop_na(x)
            dat[, x] <- scale(dat[, x])
            m <- glm(f, data = dat)
            broom::tidy(m) %>% 
                filter(term == 'scores') %>%
                mutate(term = x, conf.up = estimate + 1.96 * std.error, conf.low = estimate - 1.96 * std.error) %>%
                dplyr::select(term, estimate, conf.low, conf.up, p = p.value)
        }
    )
    write.csv(odds_ratio_table, file.path(outdir, 'plt_measures', vname, 'odds_ratio.csv'))
    or_plot <- ggplot(odds_ratio_table, aes(x = estimate, y = term, color = p < 0.05)) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.up)) +
        lims(x = c(-0.5, 0.5)) +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        theme_matt(16) +
        theme(legend.position = 'bottom') +
        labs(x = 'Beta Coefficient', y = NULL, title = NULL, color = 'p < 0.05')
    ggsave(file.path(outdir, 'plt_measures', vname, 'odds_ratio.pdf'), or_plot, width = 12, height = 8)

    # Now let's look at it continuously with the derivation samples
    odds_ratio_table <- purrr::map_df(
        vars, 
        function(x) {
            f <- as.formula(paste0(x, ' ~ scores'))
            dat <- stats_df %>% drop_na(x, labels)
            dat[, x] <- scale(dat[, x])
            m <- glm(f, data = dat)
            broom::tidy(m) %>% 
                filter(term == 'scores') %>%
                mutate(term = x, conf.up = estimate + 1.96 * std.error, conf.low = estimate - 1.96 * std.error) %>%
                dplyr::select(term, estimate, conf.low, conf.up, p = p.value)
        }
    )
    write.csv(odds_ratio_table, file.path(outdir, 'plt_measures', vname, 'odds_ratio_deriv.csv'))
    or_plot <- ggplot(odds_ratio_table, aes(x = estimate, y = term, color = p < 0.05)) +
        geom_pointrange(aes(xmin = conf.low, xmax = conf.up)) +
        lims(x = c(-0.5, 0.5)) +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        theme_matt(16) +
        theme(legend.position = 'bottom') +
        labs(x = 'Beta Coefficient', y = NULL, title = NULL, color = 'p < 0.05')
    ggsave(file.path(outdir, 'plt_measures', vname, 'odds_ratio_deriv.pdf'), or_plot, width = 12, height = 8)

    # correlations as well
    cor_table <- purrr::map_df(
    vars, 
    ~cor.test(stats_df[, 'scores'], stats_df[, .x], method = 'spearman') %>% 
        broom::tidy() %>% 
        mutate(variable = .x) %>%
        dplyr::select(variable, r = estimate, p = p.value, method, alternative)
    )
write.csv(cor_table, file.path(outdir, 'plt_measures', vname, 'correlation.csv'))
}

#======================== EVENTS BY PRESS ========================
dir.create(file.path(outdir, 'events_by_press'), showWarnings = F)
# So now let's look at the events by PRESS
censors <- grep('censor', colnames(metadata), value = T)

# make a table of the events by press
events_table_long <- metadata %>%
    dplyr::select(all_of(censors), scores, preds) %>%
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
statsTab <- events_table_long %>%
    group_by(event) %>%
    do(
        as.data.frame(stats_table(., 'censor', 'scores'))
        )
write.csv(statsTab, file.path(outdir, 'events_by_press', 'stats_table.csv'))

# barplot of the boxplot of censored and uncensored events by press
p <- events_table_long %>%
    ggplot(aes(x = event, y = scores, fill = factor(censor))) +
    geom_boxplot() +
    theme_matt(16) +
    labs(x = NULL, y = 'PRESS Score', title = NULL, fill = 'Event') +
    theme_bw(24) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'top') +
    scale_fill_manual(values = c('blue', 'red'), labels = c('No Event', 'Event'))
ggsave(file.path(outdir, 'events_by_press', 'event_boxplot.pdf'), p, width = 22, height = 8)

# So now we should make KM curves for the high_low_scores of all censors
colnames(metadata_all) <- gsub('tte_', 'time_to_', colnames(metadata_all))
censors <- grep('censor_', colnames(metadata_all), value = T)
kmDat <- metadata_all %>%
    mutate(high_low_scores = factor(ntile(scores, 2), labels = c('LOW', 'HIGH')))
group <- 'high_low_scores'
survOut <- survival_analysis(kmDat, group, censors, file.path(outdir, 'events_by_press', 'km_curves'))
with(kmDat, table(high_low_scores, censor_MACLE2))
survOut$HR_df
# okay dont love that -- this is different than the paper
colnames(kmDat)
coxph_out <- coxph(Surv(time_to_MACLE2, censor_MACLE2) ~ preds + age + sex1 + race1 + ethnicity1 + smoking1, data = kmDat)
summary(coxph_out)

events <- gsub('censor_', '', censors)
hr_table <- hazard_ratios_table(metadata_all, 'scores', censors, censor_prefix = 'censor_', time_prefix = 'time_to_', per_sd = TRUE)
write.csv(hr_table, file.path(outdir, 'events_by_press', 'hazard_ratios.csv'))
hr_table$HR_ci_upper[hr_table$HR_ci_upper > 5] <- 5
forestPlot <- ggplot(hr_table, aes(x = estimate, y = censor, color = `p.value` < 0.05)) +
    geom_pointrange(aes(xmin = HR_ci_lower, xmax = HR_ci_upper)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    theme_matt(16) +
    lims(x = c(0, 5)) +
    theme(legend.position = 'bottom') +
    labs(x = 'Hazard Ratio Per SD', y = NULL, title = NULL, color = 'p < 0.05')
ggsave(file.path(outdir, 'events_by_press', 'hazard_ratios.pdf'), forestPlot, width = 12, height = 8)


#======================== Yuhe Survival Code ========================
# So my code is not working for the survival analysis
# I'm going to just copy in Yuhe's code as I want to
# repeat what she's done
df = metadata_all
df$score_cat = ifelse(df$scores<=median(df$scores),'<=M','>M')
df$composite= ifelse(df$censor_death==1 |df$censor_MI==1| df$censor_stroke==1|df$censor_ampu==1, 1, 0)
vars = c('time_to_death','time_to_MI','time_to_stroke','time_to_ampu')
df$min_time = apply(df[,c(vars)], 1, FUN = min, na.rm = TRUE) 
df$max_time = apply(df[,c(vars)], 1, FUN = max, na.rm = TRUE) 
df$time_to_composite = ifelse(df$composite==1,df$min_time, df$max_time)

## SO BASED ON THIS CODE LOOKS LIKE YUHE IS USING MACEampu
# m = coxph(Surv(time_to_MACEampu, censor_MACEampu)~score_cat+age_surgery+sex1+ethnicity1+bmi+diabetes1+carotid_artery_disease1+prior_stroke_mini_stroke1+cli, data=df)
# m = coxph(Surv(time_to_MACLE2, censor_MACLE2)~score_cat+age_surgery+sex1+ethnicity1+bmi+diabetes1+carotid_artery_disease1+prior_stroke_mini_stroke1+cli, data=df)
m = coxph(Surv(time_to_composite, composite)~scores+age_surgery+sex1+ethnicity1+bmi+diabetes1+carotid_artery_disease1+prior_stroke_mini_stroke1+cli, data=df)
ShowRegTable(m)

# m = coxph(Surv(time_to_MACLE, censor_MACLE)~score_cat, data=df)
# m = coxph(Surv(time_to_MACLE2, censor_MACLE2)~score_cat, data=df)
# m = coxph(Surv(time_to_composite, composite)~score_cat, data=df)
# ShowRegTable(m)

# # m = coxph(Surv(time_to_composite, composite)~score_cat2, data=df)
# # ShowRegTable(m)
# # tiff('KM2.tiff', res = 200, height = 1000, width = 1200)
# # pdf('KM1.pdf')
# fit <- survfit(Surv(time_to_composite, composite)~score_cat, data=df)
# library(survminer)
# ggsurvplot(fit, data = df,
#            pval = T, pval.coord = c(0, 0.60),
#            risk.table = T,legend.title="", survscale = 'percent',
#            risk.table.height = 0.20,censor = F, legend.labs=c('Score<=Median','Score>Median'),
#            #legend.labs = c('Tertile 1','Tertile 2','Tertile 3'),
#            xlab = c('Time in days'),ylab = c('Cumulative Disease Rate'), fun = 'event', palette = 'lancet',
#            break.time.by=365,
#            tables.theme = theme_cleantable(),tables.y.text = FALSE)#+guides(color=guide_legend(nrow=2,byrow=TRUE))
# # dev.off()

#======================== Timepoint 2 ========================
dir.create(file.path(outdir, 'timepoints'), showWarnings = F)
# PRESS Scoring at Timepoint 2
genes <- pull(read.csv('data/press451_genes.csv'))
pace_path <- '../../datasets/pace/platelet_rna_filtered/dds.rds'
pace_dds <- readRDS(pace_path)
pace_counts <- normalize_counts(pace_dds, method = 'mor', log = T) %>% 
    add_missing_rows(genes) %>% 
    .[genes, ] %>%
    t()
write.csv(pace_counts, file.path(outdir, 'timepoints', 'pace_counts_press451.csv'))
dim(pace_counts)

cmd <- glue::glue(
'python code/run_press.py --data {file.path(outdir, "timepoints", "pace_counts_press451.csv")} --out {file.path(outdir, "timepoints", "pace_press_scores_timepoint2.csv")}'
)
system(cmd)

# load the scores
press_scores_timepoint2 <- read.csv(file.path(outdir, 'timepoints', 'pace_press_scores_timepoint2.csv'), row.names = 1)
timepoint_metadata <- as.data.frame(colData(pace_dds))
timepointDat <- cbind(timepoint_metadata, press_scores_timepoint2)

# get the intersect of the press scores and the timepointDat
samples <- intersect(rownames(press_scores), rownames(timepointDat))
# timepointDat[samples, c('preds', 'scores')] <- press_scores[samples, c('preds', 'scores')]

# only keep the values with duplicates in ID
timepointDat$ID <- stringr::str_replace_all(rownames(timepointDat), '\\.3', '')
timepointDat <- timepointDat[duplicated(timepointDat$ID)|duplicated(timepointDat$ID, fromLast = T),]

timeDat <- timepointDat %>% 
    dplyr::select(ID, timepoint, scores, censor_MACEampu) %>%
    pivot_wider(names_from = timepoint, values_from = scores) %>%
    mutate(
        waterfall = followup - baseline,
        consistent = (followup > 0.38) == (baseline > 0.38)
        )
timepointDat %>% dplyr::select(epi_04um_300s_n)
p <- timeDat %>%
    ggplot(aes(x = baseline, y = followup)) +
    geom_point(aes(col = consistent)) +
    stat_cor(method = 'spearman') +
    geom_hline(yintercept = 0.38, linetype = 'dashed') +
    geom_vline(xintercept = 0.38, linetype = 'dashed') +
    ggpmisc::stat_poly_eq(vjust = 2.5) +
    geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
    geom_smooth(method = 'lm', se = TRUE) +
    theme_matt(14) +
    theme(legend.position = 'top') +
    labs(x = 'Baseline', y = 'Follow Up', title = NULL)
ggsave(file.path(outdir, 'timepoints', 'press_t1_v_t2.pdf'), p, width = 4, height = 4)
with(timeDat, table(consistent))

# waterfall plot
p <- timeDat %>%
    ggplot(aes(
        x = fct_reorder(ID, waterfall),
        y = waterfall, 
        fill = factor(censor_MACEampu, labels = c('MACEAmpu', 'No MACEAmpu'))
        )) +
    geom_bar(stat = 'identity') +
    theme_matt(14) +
    labs(x = 'ID', y = 'Waterfall', title = NULL, fill = 'Censor') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position = 'top')
ggsave(file.path(outdir, 'timepoints', 'press_waterfall.pdf'), p, width = 10, height = 6)

#======================== Bleeding ========================
dir.create(file.path(outdir, 'bleeding'), showWarnings = F)
# comp bleeding
meta <- metadata
write.csv(meta, file.path(outdir, 'bleeding', 'metadata.csv'))
bleeding_comparisons <- grep('comp_bleeding', colnames(meta), value = TRUE)

# make a table of the events by press
bleeding_table_long <- meta %>%
    dplyr::select(bleeding_comparisons, scores) %>%
    pivot_longer(cols = c(bleeding_comparisons), names_to = 'event', values_to = 'censor') %>%
    group_by(event, censor)

summary_table <- bleeding_table_long %>%
    summarise(
        mean_score = mean(scores, na.rm = T),
        sd_score = sd(scores, na.rm = T),
        n = n(),
        n_censored = sum(censor == 1, na.rm = T),
        n_uncensored = sum(censor == 0, na.rm = T)
        )
write.csv(summary_table, file.path(outdir, 'bleeding', 'summary_table.csv'))

statsTab <- bleeding_table_long %>%
    group_by(event) %>%
    do(
        as.data.frame(stats_table(., 'censor', 'scores'))
        )
write.csv(statsTab, file.path(outdir, 'bleeding', 'stats_table.csv'))

# barplot of the boxplot of censored and uncensored events by press
p <- bleeding_table_long %>% drop_na() %>%
    mutate(event = gsub('comp_bleeding__', '', event)) %>%
    ggplot(aes(x = event, y = scores, fill = factor(censor))) +
    geom_boxplot() +
    geom_point(position=position_jitterdodge()) +
    theme_matt(16) +
    labs(x = NULL, y = 'PRESS Score', title = NULL, fill = 'Event') +
    theme_bw(24) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), legend.position = 'top') +
    scale_fill_manual(values = c('blue', 'red'), labels = c('No Bleed', 'Bleed')) +
    stat_compare_means(method = 't.test')
ggplot2::ggsave(file.path(outdir, 'bleeding', 'event_boxplot.pdf'), p, width = 12, height = 8)


#======================== PAD Modules ========================
# I'd like to compare PRESS to the PAD modules that I defined
# in the PAD Phenotyping paper.
module_eigengenes <- read.table('data/pace_pad_phenotyping/module_eigengenes.txt', header = T, row.names = 1)
module_samples <- intersect(rownames(module_eigengenes), rownames(metadata_all))

meta_with_modules <- cbind(metadata_all[module_samples,], module_eigengenes[module_samples,])
module_labels <- read.csv('data/pace_pad_phenotyping/module_enrichment_labels.csv', row.names = 1)
rownames(module_labels)

corrs <- purrr::map_df(
    colnames(module_eigengenes), 
    ~cor.test(meta_with_modules$scores, meta_with_modules[, .x], method = 'spearman') %>% 
        broom::tidy() %>% 
        mutate(module = gsub('ME', '', .x)) %>%
        dplyr::select(module, estimate, p.value, method, alternative)
    )
corrs_labeled <- cbind(corrs, module_labels[corrs$module,])
write.csv(corrs_labeled, file.path(outdir, 'module_correlations.csv'))
corrs_labeled$module <- factor(corrs_labeled$module, levels = rev(rownames(module_labels)))

corP <- ggplot(corrs_labeled, aes(x = estimate, y = module, color = p.value < 0.05)) +
    geom_pointrange(aes(x = estimate, y = module, xmin = 0, xmax = estimate)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    theme_matt(16) +
    theme(legend.position = 'bottom') +
    labs(x = 'Spearman Correlation', y = NULL, title = NULL, color = 'p < 0.05')
ggsave(file.path(outdir, 'module_correlations.pdf'), corP, width = 12, height = 8)
