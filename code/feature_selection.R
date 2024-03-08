###########################################################################
#
#                            feature_selection
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: feature_selection
# Output directory:
experiment <- "feature_selection"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

# So in this case we want to take the hypergeometric test results and use them to select features
# We already have only 451 that we are trying to narrow down. We can do this my subsetting to the
# highest basemeans, highest fold changes, and the hub genes.

#======================== LIBRARIES ========================#
library(tidyverse)
library(glue) # nolint

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')

#======================== CODE ========================#
# load in the stattable
stattabl_path <- './output/hyper_geneset_creation/run18_hyper60_hypo40_AGRCONTROL/full_hypertest_stattable.csv'
stats <- read.csv(stattabl_path, row.names = 1)

#======================== Feature Selection Techniques ========================#
# Basemean Selection
dir.create(file.path(outdir, 'basemean'), showWarnings = F)
head(stats)
press451 <- './output/hyper_geneset_creation/run18_hyper60_hypo40_AGRCONTROL/custom_mgc_all_GOI.txt'
press451 <- readLines(press451)

basemean_hist <- stats[press451, ] %>% 
    pivot_longer(cols = c(hypermean, nothypermean), names_to = 'mean_type', values_to = 'mean') %>%
    ggplot(aes(x = mean, fill = mean_type)) +
    scale_x_log10() +
    geom_histogram(alpha = 0.5, position = 'identity', bins = 20) +
    theme_bw()
ggsave(file.path(outdir, 'basemean', 'basemean_hist.png'), basemean_hist, width = 10, height = 10)

# check how many have a basemean in both over 100
stats[press451, ] %>% 
    filter(hypermean > 2^10 & nothypermean > 2^10)
write.csv(stats[press451, ] %>% 
    filter(hypermean > 2^10 & nothypermean > 2^10), file.path(outdir, 'basemean', 'basemean_over_1024.csv'))
selected_genes <- stats[press451, ] %>% 
    filter(hypermean > 2^10 & nothypermean > 2^10) %>% 
    rownames()
write.csv(selected_genes, file.path(outdir, 'basemean', 'selected_genes.csv'), row.names = F)

# plot the number of genes above the basemean threshold
xs <- c(2^1, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10, 2^11, 2^12)
ys <- sapply(xs, function(x) {
    stats[press451, ] %>% 
        filter(hypermean > x & nothypermean > x) %>% 
        nrow()
})

basemean_plot <- tibble(x = xs, y = ys) %>% 
    ggplot(aes(x = x, y = y)) +
    geom_point() +
    geom_line() +
    labs(x = 'Basemean Threshold', y = 'Number of Genes') +
    theme_bw(20)
ggsave(file.path(outdir, 'basemean', 'basemean_plot.png'), basemean_plot, width = 10, height = 10)

# So based on the above histogram we can see what basemeans are not a good way to select features



##### LOG2FOLD CHANGE SELECTION #####
dir.create(file.path(outdir, 'log2foldchange'), showWarnings = F)
log2foldchange_hist <- stats[press451, ] %>% 
    pivot_longer(cols = c(hyperlog2foldchange, nothyperlog2foldchange), names_to = 'log2foldchange_type', values_to = 'log2foldchange') %>%
    ggplot(aes(x = log2foldchange, fill = log2foldchange_type)) +
    geom_histogram(alpha = 0.5, position = 'identity', bins = 20) +
    theme_bw()
colnames(stats)

# load in the press genes
pace_dge <- read.csv('data/PACE_hyper_v_hypo_deseqoutput.csv', header = TRUE, row.names = 1)
colnames(pace_dge) <- paste0('pace_', colnames(pace_dge))
duke_dge <- read.csv('data/duke_validation_run3/dge_analysis/comp_group1__hyper_v_nothyper_AGCONTROL_deseqout.csv', header = TRUE, row.names = 1)
colnames(duke_dge) <- paste0('duke_', colnames(duke_dge))
comb_dge <- cbind(pace_dge[press451, ], duke_dge[press451, ])
same_dge <- comb_dge %>%
    filter(
        pace_log2FoldChange > 0 & duke_log2FoldChange > 0 |
        pace_log2FoldChange < 0 & duke_log2FoldChange < 0
        )
write.csv(same_dge, file.path(outdir, 'log2foldchange', 'same_dge.csv'))
write.csv(comb_dge, file.path(outdir, 'log2foldchange', 'comb_dge.csv'))
write.csv(rownames(same_dge), file.path(outdir, 'log2foldchange', 'selected_genes.csv'), row.names = FALSE)

p <- ggplot(same_dge, aes(x = pace_log2FoldChange, y = duke_log2FoldChange)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_smooth(method = 'lm', se = TRUE) +
    stat_cor(method = 'pearson', size = 6) +
    labs(x = 'PACE log2FoldChange', y = 'DUKE log2FoldChange') +
    geom_text_repel(aes(label = rownames(same_dge))) +
    theme_matt()
ggsave(file.path(outdir, 'log2foldchange', 'pace_duke_same_genes_lo2FC.png'), p, width = 6, height = 6)

#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
