###########################################################################
#
#                            harp_plt_paper_press
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-03-29
# Script Name: harp_plt_paper_press
# Output directory:
experiment <- "harp_plt_paper_press"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

#======================== LIBRARIES ========================
library(tidyverse)
library(glue)
library(readxl)
library(DESeq2)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')

#======================== CODE ========================
# load in the data
harp_rawcounts <- read.csv('data/harp_plt/harp_allmicontrol__rawcounttable.csv', header = TRUE, row.names = 1)
harp_metadata <- read.csv('data/harp_plt/harp_allmicontrol__metatable.csv', header = TRUE, row.names = 1)
harp_metadata_addon <- read.csv('data/harp_plt/harp-platelet-metadata-ADDON-20220823.csv', header = TRUE, row.names = 1)
harp_metadata_addon2 <- read_xlsx('data/harp_plt/Metadata HARP for WB RNAseq 3_14_2022x.xlsx')

# subset the metadata addon to the harp metadata
harp_metadata_addon <- harp_metadata_addon[rownames(harp_metadata),]

# add the addon to the metadata
harp_metadata <- harp_metadata_addon

# gsub the harp metadata names
rownames(harp_metadata) <- gsub("-", "\\.", rownames(harp_metadata))
harp_metadata_addon2$Subject.ID <- gsub("-", "\\.", harp_metadata_addon2$Subject.ID)

# subset to the press genes
press <- read.csv('data/press451_genes.csv', header = FALSE) %>% pull(.)

# add missing rows of press genes
harp_rawcounts %>% add_missing_rows(press) -> harp_rawcounts

# make the se and save it
harp_se <- make_se(harp_rawcounts, harp_metadata)
harp_dds <- DESeqDataSet(harp_se, design = ~ 1) %>% DESeq()
# save_se(harp_dds[press,], file.path(outdir, "harp_se.RData"), normalize = 'mor', log = TRUE)


# Now from here we are going to move to python for predictions.
preds <- read.csv('output/harp_predictions__run_2/harp_hyper_v_hypo_predictions.csv', header = TRUE, row.names = 1)
harp_metadata$press <- scale(preds$harp_preds)

# plot the four groups of harp
p1 <- ggplot(harp_metadata, aes(x = Angiography.Report, y = press, fill = Angiography.Report)) +
  geom_boxplot() +
  theme_matt(18) +
  theme(legend.position = 'none') +
  stat_compare_means()
p1

# change HARP.01.0135.1 to MI-CAD
harp_metadata[rownames(harp_metadata) == 'HARP.01.0135.1', 'Angiography.Report'] <- 'MI-CAD'
# change HARP.01.5049.2 to baseline
harp_metadata[rownames(harp_metadata) == 'HARP.01.5049.2', 'timepoint'] <- 'baseline'

# merge the four groups into just MI and Control
harp_metadata_full <- harp_metadata %>%
    mutate(
        ID = gsub('.[1-2]$','', rownames(harp_metadata)),
        # groupings
        MI_v_Ctrl = case_when(
            Angiography.Report == 'MI-CAD' ~ 'MI',
            Angiography.Report == 'MINOCA' ~ 'MI',
            Angiography.Report == 'Obstructive' ~ 'Control',
            Angiography.Report == 'Non-obstructive' ~ 'Control',
            TRUE ~ NA_character_
        ),
        MICAD_v_Ctrl = case_when(
            Angiography.Report == 'MI-CAD' ~ 'MI-CAD',
            Angiography.Report == 'Obstructive' ~ 'Control',
            Angiography.Report == 'Non-obstructive' ~ 'Control',
            TRUE ~ NA_character_
        ),
        MINOCA_v_Ctrl = case_when(
            Angiography.Report == 'MINOCA' ~ 'MINOCA',
            Angiography.Report == 'Obstructive' ~ 'Control',
            Angiography.Report == 'Non-obstructive' ~ 'Control',
            TRUE ~ NA_character_
        ),
        MICAD_v_MINOCA = case_when(
            Angiography.Report == 'MI-CAD' ~ 'MI-CAD',
            Angiography.Report == 'MINOCA' ~ 'MINOCA',
            TRUE ~ NA_character_
        )
    )

harp_metadata_t2 <- harp_metadata_full %>% filter(timepoint == 'followup')
harp_metadata <- harp_metadata_full %>%
    filter(timepoint == 'baseline') %>%
    mutate(
        # press tiles
        press_tile_2 = factor(ntile(press, 2), labels = c('Low', 'High')),
        press_tile_3 = factor(ntile(press, 3), labels = c('Low', 'Medium', 'High'))
    )

# add the bmi
harp_metadata$bmi <- harp_metadata_addon2 %>%
  as.data.frame() %>%
  column_to_rownames('Subject.ID') %>%
  .[rownames(harp_metadata), 'BMI']

#======================== PRESS MINOCA PAPER FIGURES ========================
# So some figures for Tessa's minoca paper
# first we need to plot:
# 1. MI v. Controls
# 2. MICAD v. Controls
# 3. MINOCA v. Controls
# 4. MICAD v. MINOCA
dir.create(file.path(outdir, 'press_by_group'), showWarnings = F)
comparing_groups <- c('MI_v_Ctrl', 'MICAD_v_Ctrl', 'MINOCA_v_Ctrl', 'MICAD_v_MINOCA', 'Angiography.Report')
for (group in comparing_groups) {
    p1 <- ggplot(drop_na(harp_metadata, group), aes(x = !!sym(group), y = press, fill = !!sym(group))) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.5) +
        labs(x = NULL, y = "PRESS Score") +
        theme_matt(18) +
        theme(legend.position = 'none') +
        scale_fill_manual(values = cmapping) +
        stat_compare_means(method = 't.test')
    ggsave(file.path(outdir, glue('press_by_group/harp_predictions_{group}.pdf')), p1)
}

dir.create(file.path(outdir, 'press_by_tile'), showWarnings = F)
tiles <- c('press_tile_2', 'press_tile_3')
cmapping <- c('MI-CAD' = 'red', 'MINOCA' = 'orange', 'Obstructive' = 'blue', 'Non-obstructive' = 'lightblue', 'MI' = 'darkred', 'Control' = 'darkblue')
for (tile in tiles) {
    for (group in comparing_groups) {
        fisher <- harp_metadata %>% 
            drop_na(group) %>% 
            dplyr::select(!!sym(tile), !!sym(group)) %>%
            table() %>%
            fisher.test()
        p1 <- harp_metadata %>%
            drop_na(group) %>%
            group_by(!!sym(tile), !!sym(group)) %>%
            summarise(n = n()) %>%
            group_by(!!sym(tile)) %>%
            mutate(pct = n / sum(n)) %>%
            ggplot(aes(x = !!sym(tile), y = pct, fill = !!sym(group))) +
            geom_bar(stat = 'identity') +
            labs(x = 'PRESS Tile', y = 'Proportion of Subjects') +
            theme_classic(18) +
            theme(legend.position = 'top') +
            scale_fill_manual(values = cmapping) +
            annotate('text', x = 1, y = 0.9, label = paste('fisher, p =', signif(fisher$p.value, 3)), size = 6)
        ggsave(file.path(outdir, glue('press_by_tile/harp_predictions_{tile}_{group}.pdf')), p1)
    }
}

dir.create(file.path(outdir, 'baseline_followup'), showWarnings = F)
for (group in comparing_groups) {
    p1 <- harp_metadata_full %>%
        drop_na(group) %>%
        group_by(ID) %>%
        filter(all(c("baseline", "followup") %in% timepoint)) %>%
        ggplot(aes(x = timepoint, y = press)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(aes(color = !!sym(group))) +
        geom_line(aes(group = ID)) +
        labs(x = 'PRESS Score', y = NULL) +
        theme_matt(18) +
        theme(legend.position = 'none') +
        stat_compare_means(method = 't.test', paired = TRUE) +
        facet_wrap(vars(!!sym(group)))
    ggsave(file.path(outdir, glue('baseline_followup/harp_predictions_{group}.pdf')), p1)
    wide_data <- harp_metadata_full %>%
        drop_na(group) %>%
        group_by(ID) %>%
        filter(all(c("baseline", "followup") %in% timepoint)) %>%
        dplyr::select(!!sym(group), press, timepoint) %>%
        pivot_wider(names_from = timepoint, values_from = press)
    p2 <- ggplot(wide_data, aes(x = baseline, y = followup)) +
        geom_point(aes(color = !!sym(group))) +
        geom_smooth(method = 'lm', se = T) +
        ggpmisc::stat_poly_eq(formula = y ~ x, vjust = 3) +
        stat_cor(method = 'pearson') +
        labs(x = 'Baseline PRESS Score', y = 'Followup PRESS Score') +
        theme_matt(18) +
        theme(legend.position = 'top') +
        scale_fill_manual(values = c('press' = 'red', 'timepoint' = 'blue')) +
        facet_wrap(vars(!!sym(group)))
    ggsave(file.path(outdir, glue('baseline_followup/harp_predictions_{group}_scatter.pdf')), p2, width = 6, height = 4)
}


#======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
