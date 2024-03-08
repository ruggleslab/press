###########################################################################
#
#                            chord_data_cleaning
#
###########################################################################
# Author: Matthew Muller
# Date: 2023-12-23
# Script Name: chord_data_cleaning
# Output directory:
experiment <- "chord_data_cleaning"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

#======================== LIBRARIES ========================#
packages <- c("lintr", "httpgd", "languageserver", "devtools", "sys", "dplyr", "tidyverse", "ggplot2", "ggpubr", "SummarizedExperiment")
pkgs <- lapply(packages, function(x) suppressMessages(require(x, character.only=T,quietly=T))) # nolint
print(packages[!sapply(pkgs, isTRUE)])

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/rna_functions.R')

#======================== CODE ========================#
# load in the data
chord_dds <- readRDS('data/chord/chord_preprocessing/dds.rds')
gene_convert <- read.csv('data/chord/gene_convert.csv', header = TRUE, row.names = 1)

# swap in the gene names
rownames(chord_dds) <- gene_convert[rownames(chord_dds), 'gene_name']

# Okay now we need to subset the chord data to the press genes
press <- read.csv('data/press451_genes.csv', header = T) %>% pull(.)

# subset the chord data
# filter to only the first time point
# chord_dds <- chord_dds[, chord_dds$Study_Timepoint == 1]

# normalize the counts
chord_dds <- DESeq(chord_dds)
chord_counts <- normalize_counts(chord_dds, method = 'log2-mor')
# chord_counts <- singscore::rankGenes(chord_dds)

# the other option is to impute the values if they are missing...
# let's do this for now
chord_counts <- add_missing_rows(chord_counts, press)

# save the chord counts
write.csv(chord_counts[press, ], file.path(outdir, "chord_press_counts.csv"))

# Fix some variables
chord_dds$Gender <- factor(chord_dds$Gender, labels = c('Male', 'Female'))

# So now we are loading in the scores and comparisons from 'chord_press_scores.ipynb'
chord_data_scoring <- read.csv('output/chord_press_scores/chord_data_with_press_scores.csv', header = TRUE, row.names = 1)
head(chord_data_scoring)

comparisons <- c('HHvLL', 'HH.HMvML.LL', 'Universal_HvL', 'UHvLL', 'HHvUL')
comparisons_clean <- c('High-high vs Low-low', 'High-high.High-mid vs Mid-low.Low-low', 'Universal high vs low', 'Universal high vs low-low', 'High-high vs Universal low') # nolint

# set the names to clean
colnames(chord_data_scoring)[colnames(chord_data_scoring) %in% comparisons] <- comparisons_clean

# stats table
# set 'nan' to NA
chord_data_scoring[chord_data_scoring == 'nan'] <- NA
chord_data_scoring$press <- factor(chord_data_scoring$press, labels = c('PRESS Low', 'PRESS High'))
chord_stats <- stats_table(chord_data_scoring, 'press', comparisons_clean, printArgs = list(showAllLevels = TRUE, printToggle = FALSE)) 
write.csv(chord_stats, file.path(outdir, "chord_press_stats.csv"))

library(magrittr)
# So not seeing much here, let's look at the breakdowns of the groupings
metadata <- as.data.frame(colData(chord_dds)) %>% filter(Study_Timepoint == 1)
rownames(metadata) <- gsub('\\.1', '', rownames(metadata))
metadata <- cbind(metadata, chord_data_scoring)

# get the epi data
epiDat <- readxl::read_xlsx('data/chord/CHORD Final Hyper Hypo Subjects copy.xlsx', sheet = 'Sheet1')
epiDat %<>% column_to_rownames('record_id')
metadata <- cbind(metadata, epiDat[rownames(metadata), ])
metadata %<>% mutate_at(vars(colnames(epiDat)), as.numeric)

#======================== Plot the scores ========================#
# plot the scores
plot_df <- metadata %>%
    dplyr::select(press_score, press, all_of(comparisons_clean)) %>%
    pivot_longer(cols = -c(press_score, press), names_to = 'comparison', values_to = 'group') %>%
    drop_na() %>%
    ggboxplot(x = 'group', y = 'press_score', fill = 'group', palette = 'jco', facet.by = 'comparison', outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    stat_compare_means(vjust = 0.75, color = 'red')
ggsave(file.path(outdir, "chord_press_scores.png"), plot = plot_df, width = 8, height = 8, dpi = 300)

#======================== Breakdowns of Groupings ========================#
# correlation plots of the epi measures to the press_score
epi_vars <- grep(colnames(epiDat), pattern = '!?lag', invert = TRUE, value = TRUE)
vars <- c("Age_at_Consent", "Gender", "Race", "BMI", "Ethnicity", "press_score")
press_stats_table <- stats_table(metadata, 'press', c(vars, comparisons_clean, epi_vars), printArgs = list(showAllLevels = TRUE, printToggle = FALSE))
write.csv(press_stats_table, file.path(outdir, "chord_press_stats_breakdown.csv"))

colnames(metadata)
for (comparison in comparisons_clean){
    vars <- c("Age_at_Consent", "Gender", "Race", "BMI", "Ethnicity", "press_score", "press")
    statTab <- stats_table(metadata, comparison, c(vars, epi_vars), printArgs = list(showAllLevels = TRUE, printToggle = FALSE))
    write.csv(statTab, file.path(outdir, glue::glue("chord_press_stats_breakdown_{comparison}.csv")))
}


plot_cors <- metadata %>%
    dplyr::select(press_score, all_of(epi_vars)) %>%
    pivot_longer(cols = -c(press_score), names_to = 'variable', values_to = 'value') %>%
    drop_na() %>%
    ggplot(aes(x = value, y = press_score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) +
    facet_wrap(~variable, scales = 'free') +
    stat_cor(method = 'spearman', color = 'red') +
    theme_bw() +
    labs(x = 'Value', y = 'PRESS Score')
ggsave(file.path(outdir, "chord_press_score_correlations.png"), plot = plot_cors, width = 9, height = 9, dpi = 300)

#======================== Plot CHORD PRESS ========================#
# make a new SE
counts <- as.data.frame(assay(chord_dds)) %>% add_missing_rows(., press)
colnames(counts) <- gsub('\\.1', '', colnames(counts))
dds <- DESeqDataSet(make_se(counts, metadata), design = ~press)

dds$HH_LL <- dds$`High-high vs Low-low`
dds$UH_UL <- dds$`Universal high vs low`

# chord_press_hm <- plot_gene_heatmap(
#     dds[press, ], 
#     "PRESS GENESET (N=451) in CHORD", 
#     c('press_score', 'HH_LL', 'UH_UL'),
#     normalize = 'vst',
#     show_column_names = FALSE,
#     show_row_names = FALSE
#     )
# pdf(file.path(outdir, "chord_press_heatmap.pdf"), width = 10, height = 10)
# print(chord_press_hm)
# dev.off()

#======================== Differential Expression ========================#
# DGE from the chord data
colnames(colData(dds)) <- gsub(' |-', '_', colnames(colData(dds)))

# set the factors correctly
conditions <- gsub(' |-', '_', comparisons_clean)
for (condition in conditions){
    dds[[condition]] <- factor(dds[[condition]], levels = c('Low', 'High'))
}

controls <- c('age_cat', 'Gender', 'Race', 'Ethnicity')
deseq_out <- deseq_analysis(dds, conditions, controls, outdir)

# load in the pace dge
pace_dge <- read.csv('data/PACE_hyper_v_hypo_deseqoutput.csv', header = TRUE, row.names = 1)
duke_dge <- read.csv('data/duke_validation_run3/dge_analysis/comp_group1__hyper_v_nothyper_AGCONTROL_deseqout.csv', header = TRUE, row.names = 1)
head(duke_dge)

# merge the two
x = as.data.frame(deseq_out[['High_high_vs_Low_low']])
colnames(x) <- paste0(colnames(x), '.x')
y = as.data.frame(pace_dge)
colnames(y) <- paste0(colnames(y), '.y')
z = as.data.frame(duke_dge)
colnames(z) <- paste0(colnames(z), '.z')
genes <- intersect(rownames(x), rownames(y))
genes <- intersect(genes, rownames(z))
dge_df <- cbind(x[press, ], y[press, ], z[press, ])
head(dge_df)

# plot the correlation
plot_df <- dge_df %>%
    ggplot(aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) +
    stat_cor(method = 'spearman', color = 'red', vjust = -0.75) +
    theme_bw() +
    labs(x = 'CHORD PRESS log2FC', y = 'PACE log2FC') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_text_repel(aes(label = rownames(dge_df[press, ])), size = 6) +
    # add labels for the number of genes in each quadrant
    annotate("text", x = 1.25, y = 2, label = paste0('Q1: ', sum(dge_df$log2FoldChange.x > 0 & dge_df$log2FoldChange.y > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = 2, label = paste0('Q2: ', sum(dge_df$log2FoldChange.x < 0 & dge_df$log2FoldChange.y > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = -2, label = paste0('Q3: ', sum(dge_df$log2FoldChange.x < 0 & dge_df$log2FoldChange.y < 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = 1.25, y = -2, label = paste0('Q4: ', sum(dge_df$log2FoldChange.x > 0 & dge_df$log2FoldChange.y < 0, na.rm = TRUE)), size = 6, color = 'red')
ggsave(file.path(outdir, "chord_press_pace_correlation.png"), plot = plot_df, width = 9, height = 9, dpi = 300)

# now do the same for the duke data
# plot the correlation
plot_df <- dge_df %>%
    ggplot(aes(x = log2FoldChange.x, y = log2FoldChange.z)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) +
    stat_cor(method = 'spearman', color = 'red', vjust = -0.75) +
    theme_bw() +
    labs(x = 'CHORD PRESS log2FC', y = 'DUKE log2FC') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_text_repel(aes(label = rownames(dge_df[press, ])), size = 6) +
    # add labels for the number of genes in each quadrant
    annotate("text", x = 1.25, y = 2, label = paste0('Q1: ', sum(dge_df$log2FoldChange.x > 0 & dge_df$log2FoldChange.z > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = 2, label = paste0('Q2: ', sum(dge_df$log2FoldChange.x < 0 & dge_df$log2FoldChange.z > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = -2, label = paste0('Q3: ', sum(dge_df$log2FoldChange.x < 0 & dge_df$log2FoldChange.z < 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = 1.25, y = -2, label = paste0('Q4: ', sum(dge_df$log2FoldChange.x > 0 & dge_df$log2FoldChange.z < 0, na.rm = TRUE)), size = 6, color = 'red')
ggsave(file.path(outdir, "chord_press_duke_correlation.png"), plot = plot_df, width = 9, height = 9, dpi = 300)

# now do the same for PACE and DUke
plot_df <- dge_df %>%
    ggplot(aes(x = log2FoldChange.y, y = log2FoldChange.z)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) +
    stat_cor(method = 'spearman', color = 'red', vjust = -0.75) +
    theme_bw() +
    labs(x = 'PACE log2FC', y = 'DUKE log2FC') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_text_repel(aes(label = rownames(dge_df[press, ])), size = 6) +
    # add labels for the number of genes in each quadrant
    annotate("text", x = 1.25, y = 2, label = paste0('Q1: ', sum(dge_df$log2FoldChange.y > 0 & dge_df$log2FoldChange.z > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = 2, label = paste0('Q2: ', sum(dge_df$log2FoldChange.y < 0 & dge_df$log2FoldChange.z > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = -2, label = paste0('Q3: ', sum(dge_df$log2FoldChange.y < 0 & dge_df$log2FoldChange.z < 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = 1.25, y = -2, label = paste0('Q4: ', sum(dge_df$log2FoldChange.y > 0 & dge_df$log2FoldChange.z < 0, na.rm = TRUE)), size = 6, color = 'red')
ggsave(file.path(outdir, "pace_duke_correlation.png"), plot = plot_df, width = 9, height = 9, dpi = 300)

# All genes now between Pace and Chord
dge_df <- cbind(x[genes, ], y[genes, ])
# plot the correlation
plot_df <- dge_df %>%
    ggplot(aes(x = log2FoldChange.x, y = log2FoldChange.y)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) +
    stat_cor(method = 'spearman', color = 'red', vjust = -0.75) +
    theme_bw() +
    labs(x = 'CHORD PRESS log2FC', y = 'PACE log2FC') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_text_repel(aes(label = rownames(dge_df[genes, ])), size = 6) +
    # add labels for the number of genes in each quadrant
    annotate("text", x = 1.25, y = 2, label = paste0('Q1: ', sum(dge_df$log2FoldChange.x > 0 & dge_df$log2FoldChange.y > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = 2, label = paste0('Q2: ', sum(dge_df$log2FoldChange.x < 0 & dge_df$log2FoldChange.y > 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = -1, y = -2, label = paste0('Q3: ', sum(dge_df$log2FoldChange.x < 0 & dge_df$log2FoldChange.y < 0, na.rm = TRUE)), size = 6, color = 'red') +
    annotate("text", x = 1.25, y = -2, label = paste0('Q4: ', sum(dge_df$log2FoldChange.x > 0 & dge_df$log2FoldChange.y < 0, na.rm = TRUE)), size = 6, color = 'red')
ggsave(file.path(outdir, "chord_press_pace_correlation_all.png"), plot = plot_df, width = 9, height = 9, dpi = 300)

#======================== Look at a subset of patients we are confident about ========================#
dir.create(file.path(outdir, 'subset'), showWarnings = FALSE)
# I went through the CHORD data…..
# I don’t recall – but please do what we first discussed
# After that – can you please compare PRESS in
# Patients with CHORD 1, 3, 13, and 21  (LOW)
# VERSUS
# CHORD 6, 8, 11, 12, 15, 18 and 19 (HIGH)
ppl <- c(
    # LOW
    'CRD001', 'CRD003', 'CRD013', 'CRD021',
    # HIGH
    'CRD006', 'CRD008', 'CRD011', 'CRD012', 'CRD015', 'CRD018', 'CRD019'
    )

tmpDat <- metadata[ppl, ]
# save the data for these people
write.csv(tmpDat, file.path(outdir, 'subset', "chord_press_scores_subset.csv"))

# plot the scores
plot_df <- tmpDat %>%
    dplyr::select(press_score, press, all_of(comparisons_clean)) %>%
    pivot_longer(cols = -c(press_score, press), names_to = 'comparison', values_to = 'group') %>%
    drop_na() %>%
    ggboxplot(x = 'group', y = 'press_score', fill = 'group', palette = 'jco', facet.by = 'comparison', outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.5) +
    stat_compare_means(vjust = 0.75, color = 'red')
ggsave(file.path(outdir, 'subset', "chord_press_scores_subset.png"), plot = plot_df, width = 8, height = 8, dpi = 300)

# correlation plots of the epi measures to the press_score
plot_cors <- tmpDat %>%
    dplyr::select(press_score, all_of(epi_vars)) %>%
    pivot_longer(cols = -c(press_score), names_to = 'variable', values_to = 'value') %>%
    drop_na() %>%
    ggplot(aes(x = value, y = press_score)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = 'lm', se = TRUE) +
    facet_wrap(~variable, scales = 'free') +
    stat_cor(method = 'spearman', color = 'red') +
    theme_bw() +
    labs(x = 'Value', y = 'PRESS Score')
ggsave(file.path(outdir, 'subset', "chord_press_score_correlations_subset.png"), plot = plot_cors, width = 9, height = 9, dpi = 300)

#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))


