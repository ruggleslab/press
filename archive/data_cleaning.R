###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-02-15
# 
# Script Name: Data Wrangling Script in R
# 
# Notes:
# This is a general script for cohort variation analysis in R for the Duke cohort


###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "tidyverse",
  "ggplot2",
  "devtools",
  "BiocManager",
  "SummarizedExperiment",
  "gridExtra",
  "reshape2",
  "plyr",
  "dplyr",
  "devtools"
)

for (pkg in packages) {
  paste0(pkg)[[1]]
  library(pkg, character.only = T, quietly = T)
}

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url("https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R")
source_url("https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/deseq_helper_functions.R")

# Output directory
experiment <- 'data_cleaning'
outdir <- file.path('output', experiment)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)



###########################################################################
#
#                                 CODE
#
###########################################################################
# Add code here
#

duke_metadata = read.csv('data/GSE158765/MetaData_samples.csv') %>%
  dplyr::select(#'sample_name',
         'characteristic__subject_id',
         'characteristic__treatment_drug_label',
         #'characteristic__passed_sample_QC',
         #'characteristic__treatment_drug',
         'characteristic__visit_number'
         )


# crosstab <- crosstable(as.data.frame(sapply(duke_metadata, as.character)), c(characteristic__visit_number), 
#                        by=c(characteristic__treatment_drug_label),
#                        label=F, num_digits=2) %>% as_flextable(compact=T)

# duke_metadata_tosh <- read.csv('data/duke_validation_run3/dukemetatable_sel.csv')

# duke_metadata_tosh$cohort <- mapvalues(duke_metadata_tosh$cohort, 
#                                 from=unique(duke_metadata_tosh$cohort),
#                                 to=c(unique(duke_metadata_tosh$characteristic__treatment_drug_label)[-1],"washout")
#                                 )

# prep the full pace count data
pace_counts <- read.table('data/pace/rnaseq/pace_platelet_quant.featurecounts.counts.unstr.txt', header=T, row.names=1, sep='\t')
pace_counts <- read.table('data/plt_filtrawcounttab.txt', header=T, row.names=1, sep='\t')
colnames(pace_counts) <- gsub('_', '', colnames(pace_counts))
pace_metadata <- read.csv('data/pace/pace_metadata_2023-03-31.csv', row.names=1)

# load in press geneset
press_geneset <- read.csv('data/clean/press_genes.csv', header=F, stringsAsFactors = F) %>% pull()

# read in the possible outliers table
possible_outliers <- read.csv('docs/possible_outlier_tab_2.csv', stringsAsFactors = F)

# remove them from the counts
pace_counts <- pace_counts[, !(colnames(pace_counts) %in% possible_outliers$samples)]

# make the se
pace_se <- make_se(pace_counts, pace_metadata)
pace_metadata[gsub('_', '', possible_outliers$samples), ] %>% select(epi_01um_300s_n)

# make dds
pace_dds <- DESeqDataSet(pace_se, design = ~ 1)

# run preprocessing
pace_dds <- pace_dds %>% DESeq()

# get the normalized counts
pace_norm <- normalize_counts(pace_dds)

# get the normalized counts of the press genes
pace_norm_press <- pace_norm[press_geneset,]

# save the normalized counts
write.csv(pace_norm_press, 'data/clean/pace_norm_all.csv')
dim(pace_norm_press)
