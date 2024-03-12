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
library(tidyverse)
library(ggpubr)
library(ggrepel)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')

#======================== DATA ========================#
# load in the press genes
pace_dge <- read.csv('data/PACE_hyper_v_hypo_deseqoutput.csv', header = TRUE, row.names = 1)
duke_dge <- read.csv('data/duke_validation_run3/dge_analysis/comp_group1__hyper_v_nothyper_AGCONTROL_deseqout.csv', header = TRUE, row.names = 1)

# comparisons of interest
press <- pull(read.csv('data/press451_genes.csv', header = TRUE))
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

#======================== CHORD MISSING GENES ========================#
chord_missing_genes <- rownames(chord_counts)[rowSums(chord_counts) == 0]

# table of the lo2FC of the missing genes from PACE and DUKE
missing_genes_lo2FC <- data.frame(
  row.names = chord_missing_genes,
  pace = pace_dge[chord_missing_genes, 'log2FoldChange'],
  duke = duke_dge[chord_missing_genes, 'log2FoldChange']
)

p <- ggplot(missing_genes_lo2FC, aes(x = pace, y = duke)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'PACE log2FoldChange', y = 'DUKE log2FoldChange') +
  geom_text_repel(aes(label = rownames(missing_genes_lo2FC))) +
  theme_matt()
ggsave(file.path(outdir, 'pace_duke_missing_genes_lo2FC.png'), p, width = 6, height = 6)

chord_present_genes <- rownames(chord_counts)[rowSums(chord_counts) > 0]
chord_present_genes_lo2FC <- data.frame(
  row.names = chord_present_genes,
  pace = pace_dge[chord_present_genes, 'log2FoldChange'],
  duke = duke_dge[chord_present_genes, 'log2FoldChange']
)

p <- ggplot(chord_present_genes_lo2FC, aes(x = pace, y = duke)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'PACE log2FoldChange', y = 'DUKE log2FoldChange') +
  geom_text_repel(aes(label = rownames(chord_present_genes_lo2FC))) +
  theme_matt()
ggsave(file.path(outdir, 'pace_duke_present_genes_lo2FC.png'), p, width = 6, height = 6)



# now do the same for the base mean rather than log2FoldChange
missing_genes_base_mean <- data.frame(
  row.names = chord_missing_genes,
  pace = pace_dge[chord_missing_genes, 'baseMean'],
  duke = duke_dge[chord_missing_genes, 'baseMean']
)

p <- ggplot(missing_genes_base_mean, aes(x = pace, y = duke)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'PACE baseMean', y = 'DUKE baseMean') +
  geom_text_repel(aes(label = rownames(missing_genes_base_mean))) +
  theme_matt()

ggsave(file.path(outdir, 'pace_duke_missing_genes_baseMean.png'), p, width = 6, height = 6)

present_genes_base_mean <- data.frame(
  row.names = chord_present_genes,
  pace = pace_dge[chord_present_genes, 'baseMean'],
  duke = duke_dge[chord_present_genes, 'baseMean']
)

p <- ggplot(present_genes_base_mean, aes(x = pace, y = duke)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'PACE baseMean', y = 'DUKE baseMean') +
  geom_text_repel(aes(label = rownames(present_genes_base_mean))) +
  theme_matt()
ggsave(file.path(outdir, 'pace_duke_present_genes_baseMean.png'), p, width = 6, height = 6)


#======================== CHORD Timepoint Analysis ========================#
dir.create(file.path(outdir, 'chord_timepoint_analysis'), showWarnings = F)
# So would be helpful to look at the timepoints for chord
chord_t1 <- chord_scores
chord_dds <- readRDS('data/chord/chord_preprocessing/dds.rds')
gene_convert <- read.csv('data/chord/gene_convert.csv', row.names = 1)

chord_dds <- DESeq2::estimateSizeFactors(chord_dds)
chord_t2_dds <- chord_dds[, chord_dds$Study_Timepoint == 2]
rownames(chord_t2_dds) <- gene_convert[rownames(chord_t2_dds), 'gene_name']
chord_counts_t2 <- normalize_counts(chord_t2_dds, 'mor', log2 = T) %>%
  add_missing_rows(press) %>%
  .[press,] %>%
  t()
write.csv(chord_counts_t2, file.path(outdir, 'chord_timepoint_analysis', 'chord_counts_t2.csv'))

inpath <- 'output/gene_qc_investigation/chord_timepoint_analysis/chord_counts_t2.csv'
outpath <- 'output/gene_qc_investigation/chord_timepoint_analysis/chord_t2_press_scores.csv'
cmd <- glue('python3 code/run_press.py --data {inpath} --out {outpath}')
system(cmd)

chord_t2 <- read.csv(outpath, row.names = 1)

rownames(chord_t2) <- gsub('\\.2', '', rownames(chord_t2))
sizeFact_t1 <- chord_dds$sizeFactor[chord_dds$Study_Timepoint == 1]
sizeFact_t2 <- chord_dds$sizeFactor[chord_dds$Study_Timepoint == 2]
chord_t1$sizeFact_t1 <- sizeFact_t1
chord_t2$sizeFact_t2 <- sizeFact_t2
timepointDat <- cbind(chord_t1, t2_scores = chord_t2$scores) %>%
  rename('t1_scores' = 'press_score')

colnames(timepointDat)
tPlot <- ggplot(timepointDat, aes(x = t1_scores, y = t2_scores)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  ggpmisc::stat_poly_eq(size = 6, vjust = 3.5) +
  labs(x = 'Timepoint 1 Press Score', y = 'Timepoint 2 Press Score') +
  geom_text_repel(aes(label = rownames(timepointDat))) +
  theme_matt()
ggsave(file.path(outdir, 'chord_timepoint_analysis', 'chord_timepoint_timepoints.png'), tPlot, width = 6, height = 6)

tPlot <- ggplot(timepointDat, aes(x = t1_scores, y = t2_scores, color = factor(press))) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) +
  ggpmisc::stat_poly_eq(size = 6, vjust = 3.5) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'Timepoint 1 Press Score', y = 'Timepoint 2 Press Score') +
  geom_text_repel(aes(label = rownames(timepointDat))) +
  theme_matt() +
  theme(legend.position = 'top')
ggsave(file.path(outdir, 'chord_timepoint_analysis', 'chord_timepoint_analysis.png'), tPlot, width = 6, height = 6)

t1Lib <- ggplot(timepointDat, aes(x = libSize_t1, y = t1_scores)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'Timepoint 1 Size Factor', y = 'Timepoint 1 Press Score') +
  geom_text_repel(aes(label = rownames(timepointDat))) +
  theme_matt()
ggsave(file.path(outdir, 'chord_timepoint_analysis', 'chord_timepoint_t1_libsize.png'), t1Lib, width = 6, height = 6)

t2Lib <- ggplot(timepointDat, aes(x = sizeFact_t2, y = t2_scores)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'Timepoint 2 Size Factor', y = 'Timepoint 2 Press Score') +
  geom_text_repel(aes(label = rownames(timepointDat))) +
  theme_matt()
ggsave(file.path(outdir, 'chord_timepoint_analysis', 'chord_timepoint_t2_libsize.png'), t2Lib, width = 6, height = 6)

sizeFacts <- ggplot(timepointDat, aes(x = sizeFact_t2, y = sizeFact_t1)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE) +
  stat_cor(method = 'spearman', size = 6) +
  labs(x = 'Timepoint 2 Size Factor', y = 'Timepoint 1 Size Factor') +
  geom_text_repel(aes(label = rownames(timepointDat))) +
  theme_matt()
ggsave(file.path(outdir, 'chord_timepoint_analysis', 'chord_timepoint_libsize.png'), sizeFacts, width = 6, height = 6)



#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
