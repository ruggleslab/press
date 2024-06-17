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
outdir <- file.path("output", paste0(experiment))
dir.create(outdir, showWarnings = F)

#======================== LIBRARIES ========================
library(tidyverse)
library(glue)
library(readxl)
library(DESeq2)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R")

#======================== CODE ========================
# load in the data
harp_rawcounts <- read.csv("data/harp_plt/harp_allmicontrol__rawcounttable.csv", header = TRUE, row.names = 1)
harp_metadata <- read.csv("data/harp_plt/harp_allmicontrol__metatable.csv", header = TRUE, row.names = 1)
harp_metadata_addon <- read.csv("data/harp_plt/harp-platelet-metadata-ADDON-20220823.csv", header = TRUE, row.names = 1)
harp_metadata_addon2 <- read_xlsx("data/harp_plt/Metadata HARP for WB RNAseq 3_14_2022x.xlsx")

# harp_rawcounts <-  harp_rawcounts %>%
#     drop_na(gene_name) %>%
#     dplyr::select(contains("HARP"))

# subset the metadata addon to the harp metadata
harp_metadata_addon <- harp_metadata_addon[rownames(harp_metadata), ]

# add the addon to the metadata
harp_metadata <- harp_metadata_addon

# gsub the harp metadata names
rownames(harp_metadata) <- gsub("-", "\\.", rownames(harp_metadata))
harp_metadata_addon2$Subject.ID <- gsub("-", "\\.", harp_metadata_addon2$Subject.ID)

# subset to the press genes
press <- read.csv("data/press451_genes.csv", header = TRUE) %>% pull(.)

# add missing rows of press genes
harp_rawcounts <- harp_rawcounts %>% add_missing_rows(press)

# make the se and save it
harp_se <- make_se(harp_rawcounts, harp_metadata)
harp_dds <- DESeqDataSet(harp_se, design = ~ 1) %>% estimateSizeFactors()
# save_se(harp_dds[press,], file.path(outdir, "harp_se.RData"), normalize = "mor", log = TRUE)


# Now from here we are going to move to python for predictions.
preds <- read.csv("output/harp_predictions__run_2/harp_hyper_v_hypo_predictions.csv", header = TRUE, row.names = 1)
harp_metadata$press <- scale(preds$harp_preds)

# plot the four groups of harp
p1 <- ggplot(harp_metadata, aes(x = Angiography.Report, y = press, fill = Angiography.Report)) +
  geom_boxplot() +
  theme_matt(18) +
  theme(legend.position = "none") +
  stat_compare_means()
p1

# change HARP.01.0135.1 to MI-CAD
harp_metadata[rownames(harp_metadata) == "HARP.01.0135.1", "Angiography.Report"] <- "MINOCA" # NO UPDATE THIS IS WRONG # UPDATE THIS NOW MINOCA AS OF 2024-04-25!

# change HARP.01.5049.2 to baseline
harp_metadata[rownames(harp_metadata) == "HARP.01.5049.2", "timepoint"] <- "baseline"

# exclusions
exclude <- c("HARP.01.0060.1", "HARP.01.0047.1", "HARP.01.0076.1", "HARP.01.0104")
harp_metadata <- harp_metadata[!rownames(harp_metadata) %in% exclude, ]


# merge the four groups into just MI and Control
harp_metadata_full <- harp_metadata %>%
    mutate(
        ID = gsub(".[1-2]$","", rownames(harp_metadata)),
        # groupings
        MI_v_Ctrl = case_when(
            Angiography.Report == "MI-CAD" ~ "MI",
            Angiography.Report == "MINOCA" ~ "MI",
            Angiography.Report == "Obstructive" ~ "Control",
            Angiography.Report == "Non-obstructive" ~ "Control",
            TRUE ~ NA_character_
        ),
        MICAD_v_Ctrl = case_when(
            Angiography.Report == "MI-CAD" ~ "MI-CAD",
            Angiography.Report == "Obstructive" ~ "Control",
            Angiography.Report == "Non-obstructive" ~ "Control",
            TRUE ~ NA_character_
        ),
        MINOCA_v_Ctrl = case_when(
            Angiography.Report == "MINOCA" ~ "MINOCA",
            Angiography.Report == "Obstructive" ~ "Control",
            Angiography.Report == "Non-obstructive" ~ "Control",
            TRUE ~ NA_character_
        ),
        MICAD_v_MINOCA = case_when(
            Angiography.Report == "MI-CAD" ~ "MI-CAD",
            Angiography.Report == "MINOCA" ~ "MINOCA",
            TRUE ~ NA_character_
        ),
        MICAD_v_MINOCA_v_Ctrl = case_when(
          Angiography.Report == "MI-CAD" ~ "MI-CAD",
          Angiography.Report == "MINOCA" ~ "MINOCA",
          Angiography.Report == "Obstructive" ~ "Control",
          Angiography.Report == "Non-obstructive" ~ "Control",
          TRUE ~ NA_character_
        )
    )

harp_metadata_t2 <- harp_metadata_full %>% filter(timepoint == "followup")
harp_metadata <- harp_metadata_full %>%
    filter(timepoint == "baseline") %>%
    mutate(
        # press tiles
        press_tile_2 = factor(ntile(press, 2), labels = c("Low", "High")),
        press_tile_3 = factor(ntile(press, 3), labels = c("Low", "Medium", "High"))
    )

# add the bmi
harp_metadata$bmi <- harp_metadata_addon2 %>%
  as.data.frame() %>%
  column_to_rownames("Subject.ID") %>%
  .[rownames(harp_metadata), "BMI"]

#======================== PRESS MINOCA PAPER FIGURES ========================
# So some figures for Tessa"s minoca paper
# first we need to plot:
# 1. MI v. Controls
# 2. MICAD v. Controls
# 3. MINOCA v. Controls
# 4. MICAD v. MINOCA

# let's use the same data as florencia did in the paper to make sure we are on the same page
load("data/harp_plt/HARP_platelet_raw_counts_April2024.rdata", verbose = TRUE)
metadata <- read.csv("data/harp_plt/metadata_fixed_April2024.csv", header = TRUE, row.names = 1)
rownames(metadata) <- metadata$SubjectID

# prep the data for the press scores
rownames(gene_convert) <- gene_convert$ensembl_gene_id
raw_counts <- raw_counts %>%
    rownames_to_column("gene_name") %>%
    mutate(gene_name = gsub("\\..*", "", gene_name)) %>%
    mutate(gene_name = gene_convert[gene_name, "external_gene_name"]) %>%
    drop_na(gene_name) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    filter(gene_name != "") %>%
    column_to_rownames("gene_name") %>%
    add_missing_rows(press)
harp_se <- make_se(raw_counts, metadata)
harp_dds <- DESeqDataSet(harp_se, design = ~ 1) %>% estimateSizeFactors()

harp_press_counts <- normalize_counts(harp_dds, "mor", log = TRUE)[press, ]
rownames(harp_press_counts) <- press
harp_in <- file.path(outdir, "harp_press_counts.csv")
write.csv(t(harp_press_counts), harp_in)
harp_out <- file.path(outdir, "harp_press_predictions.csv")
system(glue("python code/run_press.py --data {harp_in} --out {harp_out}"))
harp_press_preds <- read.csv(harp_out, header = TRUE, row.names = 1)
harp_press_preds$scores <- scale(harp_press_preds$preds)
meta_with_press <- metadata %>%
    mutate(
        press = harp_press_preds$score,
        MI_v_Ctrl = case_when(
            Angiography.Report == "MI-CAD" ~ "MI",
            Angiography.Report == "MINOCA" ~ "MI",
            Angiography.Report == "Obstructive" ~ "Control",
            Angiography.Report == "Non-obstructive" ~ "Control",
            TRUE ~ NA_character_
        ),
        MICAD_v_Ctrl = case_when(
            Angiography.Report == "MI-CAD" ~ "MI-CAD",
            Angiography.Report == "Obstructive" ~ "Control",
            Angiography.Report == "Non-obstructive" ~ "Control",
            TRUE ~ NA_character_
        ),
        MINOCA_v_Ctrl = case_when(
            Angiography.Report == "MINOCA" ~ "MINOCA",
            Angiography.Report == "Obstructive" ~ "Control",
            Angiography.Report == "Non-obstructive" ~ "Control",
            TRUE ~ NA_character_
        ),
        MICAD_v_MINOCA = case_when(
            Angiography.Report == "MI-CAD" ~ "MI-CAD",
            Angiography.Report == "MINOCA" ~ "MINOCA",
            TRUE ~ NA_character_
        ),
        MICAD_v_MINOCA_v_Ctrl = case_when(
          Angiography.Report == "MI-CAD" ~ "MI-CAD",
          Angiography.Report == "MINOCA" ~ "MINOCA",
          Angiography.Report == "Obstructive" ~ "Control",
          Angiography.Report == "Non-obstructive" ~ "Control",
          TRUE ~ NA_character_
        ),
        press_tile_2 = factor(ntile(press, 2), labels = c("Low", "High")),
        press_tile_3 = factor(ntile(press, 3), labels = c("Low", "Medium", "High"))
    )

# I fucked up the old model so let's use the old scores
rownames(preds) <- gsub("\\.", "-", rownames(preds))
samples <- intersect(rownames(preds), rownames(meta_with_press))
meta_with_press[samples, 'press'] <- scale(preds[samples, "harp_preds"])

cmapping <- c("MI-CAD" = "red", "MINOCA" = "orange", "Obstructive" = "blue", "Non-obstructive" = "lightblue", "MI" = "darkred", "Control" = "darkblue")
comparing_groups <- c("MI_v_Ctrl", "MICAD_v_Ctrl", "MINOCA_v_Ctrl", "MICAD_v_MINOCA", "MICAD_v_MINOCA_v_Ctrl", "Angiography.Report")
tiles <- c("press_tile_2", "press_tile_3")
dat_groups <- meta_with_press %>%
    filter(timepoint == "baseline") %>%
    dplyr::select(any_of(c(comparing_groups, tiles)), press)
dat_groups$press <- signif(dat_groups$press)
write.csv(dat_groups, file.path(outdir, "harp_press_groups.csv"))

dir.create(file.path(outdir, "press_by_group"), showWarnings = FALSE)
for (group in comparing_groups) {
    comps <- pairwise_combos(na.omit(dat_groups[[group]]))
    p1 <- ggplot(drop_na(dat_groups, group), aes(x = !!sym(group), y = press, fill = !!sym(group))) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, alpha = 0.5) +
        labs(x = NULL, y = "PRESS Score") +
        theme_matt(18) +
        theme(legend.position = "none") +
        scale_fill_manual(values = cmapping) +
        stat_compare_means(method = "t.test", comparisons = comps)
    ggsave(file.path(outdir, glue("press_by_group/harp_predictions_{group}.pdf")), p1)
}

dir.create(file.path(outdir, "press_by_tile"), showWarnings = FALSE)
for (tile in tiles) {
    for (group in comparing_groups) {
        fisher <- dat_groups %>%
            drop_na(group) %>%
            dplyr::select(!!sym(tile), !!sym(group)) %>%
            table() %>%
            fisher.test()
        p1 <- dat_groups %>%
            drop_na(group) %>%
            group_by(!!sym(tile), !!sym(group)) %>%
            summarise(n = n()) %>%
            group_by(!!sym(tile)) %>%
            mutate(pct = n / sum(n)) %>%
            ggplot(aes(x = !!sym(tile), y = pct, fill = !!sym(group))) +
            geom_bar(stat = "identity") +
            labs(x = "PRESS Tile", y = "Proportion of Subjects") +
            theme_classic(18) +
            theme(legend.position = "top") +
            scale_fill_manual(values = cmapping) +
            annotate("text", x = 1, y = 0.9, label = paste("fisher, p =", signif(fisher$p.value, 3)), size = 6)
        ggsave(file.path(outdir, glue("press_by_tile/harp_predictions_{tile}_{group}.pdf")), p1)
    }
}

# dir.create(file.path(outdir, "baseline_followup"), showWarnings = F)
# for (group in comparing_groups) {
#     p1 <- harp_metadata_full %>%
#         drop_na(group) %>%
#         group_by(ID) %>%
#         filter(all(c("baseline", "followup") %in% timepoint)) %>%
#         ggplot(aes(x = timepoint, y = press)) +
#         geom_boxplot(outlier.shape = NA) +
#         geom_point(aes(color = !!sym(group))) +
#         geom_line(aes(group = ID)) +
#         labs(x = "PRESS Score", y = NULL) +
#         theme_matt(18) +
#         theme(legend.position = "none") +
#         stat_compare_means(method = "t.test", paired = TRUE) +
#         facet_wrap(vars(!!sym(group)))
#     ggsave(file.path(outdir, glue("baseline_followup/harp_predictions_{group}.pdf")), p1)
#     wide_data <- harp_metadata_full %>%
#         drop_na(group) %>%
#         group_by(ID) %>%
#         filter(all(c("baseline", "followup") %in% timepoint)) %>%
#         dplyr::select(!!sym(group), press, timepoint) %>%
#         pivot_wider(names_from = timepoint, values_from = press)
#     p2 <- ggplot(wide_data, aes(x = baseline, y = followup)) +
#         geom_point(aes(color = !!sym(group))) +
#         geom_smooth(method = "lm", se = T) +
#         ggpmisc::stat_poly_eq(formula = y ~ x, vjust = 3) +
#         stat_cor(method = "pearson") +
#         labs(x = "Baseline PRESS Score", y = "Followup PRESS Score") +
#         theme_matt(18) +
#         theme(legend.position = "top") +
#         scale_fill_manual(values = c("press" = "red", "timepoint" = "blue")) +
#         facet_wrap(vars(!!sym(group)))
#     ggsave(file.path(outdir, glue("baseline_followup/harp_predictions_{group}_scatter.pdf")), p2, width = 6, height = 4)
# }


#======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
