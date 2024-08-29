###########################################################################
#
#                            covid_gordon_figure
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-07-29
# Script Name: covid_gordon_figure
# Output directory:
experiment <- "covid_gordon_figure"
outdir <- file.path("output", experiment)
dir.create(outdir, showWarnings = FALSE)

#======================== LIBRARIES ========================
library(tidyverse)
library(glue)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R")
source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/enrichment_functions.R")

#======================== CODE ========================
# So we are going to just make a figure showing that PRESS works
# in COVID for JBs presentation at the Gordon conference
# this can be done with the data we have already
# we will want to make a 3 panels --
# 1. the PRESS score for our COVID data
# 2. the PRESS enrichment GSEA for our COVID data
# 3. the PRESS enrichment GSEA for the COVID Utah Data

press <- pull(read.csv('data/press451_genes.csv', header = T))

#======================== Panel 1 ========================
dir.create(file.path(outdir, "covid_scores"), showWarnings = FALSE)
# press for our COVID data
counts_raw <- read.csv("data/20240729_covid_figure/covid_covidcontrol/intermediate_files/covid_covidcontrol__rawcounttable.csv", row.names = 1)
meta <- read.csv("data/20240729_covid_figure/covid_covidcontrol/intermediate_files/covid_covidcontrol__metatable.csv", row.names = 1)
rownames(meta) <- gsub("-", "\\.", rownames(meta))

missing_genes <- setdiff(press, rownames(counts_raw))
counts_raw_clean <- counts_raw %>% 
    add_missing_rows(press) %>%
    filter(rownames(.) %in% press)

se <- make_se(counts_raw_clean, meta)
dds <- DESeqDataSet(se, design = ~ 1)
dds <- estimateSizeFactors(dds)
counts_logmor <- t(normalize_counts(dds, "mor", log = TRUE))
write.csv(counts_logmor, file.path(outdir, "covid_scores", "press_counts.csv"))

dat <- glue("{outdir}/covid_scores/press_counts.csv")
out <- glue("{outdir}/covid_scores/press_scores.csv")
system(glue("python3 code/run_press.py --data {dat} --out {out}"))

scores <- read.csv(file.path(outdir, "covid_scores", "press_scores.csv"), row.names = 1)
dat <- cbind(meta, scores)

covid_press_plot <- ggplot(dat, aes(x = factor(comp_covid__covid_v_control), y = scores, fill = factor(comp_covid__covid_v_control))) +
    geom_boxplot() +
    theme_minimal() +
    labs(x = "COVID vs Control", y = "PRESS Score") +
    geom_hline(yintercept = 0.38, linetype = "dashed") +
    annotate("text", x = 0.8, y = 0.38, label = "Hyperreactive", vjust = -1, color = "red", size = 6) +
    annotate("text", x = 0.8, y = 0.38, label = "Normoreactive", vjust = 1, color = "blue", size = 6) +
    theme(legend.position = "none")
ggsave(file.path(outdir, "covid_scores", "covid_press_plot.pdf"), covid_press_plot, width = 5, height = 5)

#======================== PRESS GSEA ========================
# em2 <- GSEA(geneList, TERM2GENE = C3_t2g)
# we can do this for PRESS Up and PRESS Down
press_up <- read.table("output/hyper_geneset_creation/run18_hyper60_hypo40_AGRCONTROL/custom_mgc_hyper_up.txt")
press_down <- read.table("output/hyper_geneset_creation/run18_hyper60_hypo40_AGRCONTROL/custom_mgc_hyper_down.txt")

press_up$direction <- "PRESS Up"
press_down$direction <- "PRESS Down"

press_t2g <- rbind(press_up, press_down)
colnames(press_t2g) <- c("gene", "term")
press_t2g <- press_t2g %>% dplyr::select(term, gene)

## COVID data preranked
utah_covid <- read.csv("data/20240729_covid_figure/utahcovid_covidcontrol/utahcovid_covidcontrol_deseq.csv", row.names = 1)
utah_covid_prerank <- get_fc_list(utah_covid, "log2FoldChange")

nyu_covid <- read.csv("data/20240729_covid_figure/covid_covidcontrol/covid_covidcontrol_wcontrol_deseq.csv", row.names = 1)
nyu_covid_prerank <- get_fc_list(nyu_covid, "log2FoldChange")

# get the gsea
nyu_gsea <- GSEA(nyu_covid_prerank, TERM2GENE = press_t2g)
utah_gsea <- GSEA(utah_covid_prerank, TERM2GENE = press_t2g)
covid_gsea <- rbind(
    as.data.frame(nyu_gsea) %>% mutate(dataset = "NYU"),
    as.data.frame(utah_gsea) %>% mutate(dataset = "UTAH")
)

# let's start with some ridge plots
nyu_ridge <- ridgeplot(nyu_gsea)
utah_ridge <- ridgeplot(utah_gsea)
ridges <- cowplot::plot_grid(nyu_ridge, utah_ridge, nrow = 2, labels = c("NYU", "UTAH"))
ggsave(file.path(outdir, "covid_scores", "ridge_plots.pdf"), ridges, width = 12, height = 5)

# now some walkplots
nyu_walk <- gseaplot2(nyu_gsea, 1:2)
ggplot2::ggsave(file.path(outdir, "covid_scores", "nyu_walk.pdf"), nyu_walk, width = 10, height = 5)
utah_walk <- gseaplot2(utah_gsea, 1:2)
ggplot2::ggsave(file.path(outdir, "covid_scores", "utah_walk.pdf"), utah_walk, width = 10, height = 5)

#======================== Dotplot ========================
covid_gsea[, 1:8]
dotplot <- ggplot(covid_gsea, aes(x = ID, y = dataset, size = -log(p.adjust), color = factor(sign(NES)))) +
    geom_point() +
    scale_size_area(max_size = 12) +
    scale_color_manual(values = c("blue", "red")) +
    labs(x = NULL, y = NULL, size = "-log(p.adjust)", color = "NES")
dotplot
ggplot2::ggsave(file.path(outdir, "covid_scores", "dotplot.pdf"), dotplot, width = 8, height = 4)

#======================== END ========================
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
