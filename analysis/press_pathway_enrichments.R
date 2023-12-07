###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-11
# 
# Script Name: PRESS pathway enrichments
# 
# Notes:
# This file will recreate press pathway enrichments within the pace data and
# duke as well. The goal is to create a pathway ensemble signature that can be
# used to predict duke that is built up from press


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "BiocManager",
  "devtools",
  "SummarizedExperiment",
  "clusterProfiler",
  "org.Hs.eg.db",
  "enrichplot",
  "AnnotationDbi"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T)
}

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/converting_functions.R')

# Output directory
experiment <- 'press_pathway_enrichment'
outdir <- file.path('output', paste0(experiment, '__',format(Sys.time(),"%F"), '/'))
dir.create(outdir)

###########################################################################
#
#                                 CODE
#
###########################################################################


load('data/summarized_experiments/se_pace_pltRNA.rdata')
load('data/summarized_experiments/se_duke_pltRNA.rdata')

# press without 
press <- read.csv('data/clean/press_genes_no_transcripts.csv', header = F)$V1 %>% gsub("-", ".", .)





###########################################################################
#
#                           GO Classification
#
###########################################################################

press_MFs <- groupGO(gene     = press,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "MF",
                     keyType  = "SYMBOL",
                     readable = TRUE
                     )

press_CCs <- groupGO(gene     = press,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "CC",
                     keyType  = "SYMBOL",
                     readable = TRUE
                     )

press_BPs <- groupGO(gene     = press,
                     OrgDb    = org.Hs.eg.db,
                     ont      = "BP",
                     keyType  = "SYMBOL",
                     readable = TRUE
                     )



# Testing around with enrichment
ego <- enrichGO(
  gene      = se_pace_plt,
  OrgDb     = org.Hs.eg.db,
  keyType   = "SYMBOL",
  ont       = "ALL",
  readable  = TRUE
  )

head(ego)







