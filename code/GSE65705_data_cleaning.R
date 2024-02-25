###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-04-04
# 
# Script Name: GSE65705 Data Cleaning
# 
# Notes:
# This file is general QC and datacleaning for GSE65705. This is a study focused on RNA-seq
# data to compare DGE in non-MI versus MI patients. For our purposes, we are going to
# process the data for use within an SKlearn model pipeline.

  # Output directory
experiment <- 'GSE65705_hyper_preds'
outdir <- file.path('output', paste0(experiment, '/'))
dir.create(outdir, showWarnings = F)

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
  "readxl",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "DESeq2",
  "plyr"
)
junk <- sapply(packages, function(x){library(x, character.only = TRUE, quietly = TRUE, verbose = FALSE)})


# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/scripts/main/Rtools/general_functions.R')


collect_counts_from_folder <- function(folder_path, reads_column="counts", id_var="Gene", sorted=T) {
  require(readr)
  require(dplyr)
  # Get the list of files in the folder
  file_list <- list.files(path = folder_path, full.names = TRUE)
  
  # Read the first file to create the dataframe
  initial_df <- reader(file_list[1], quiet = T) 
  combined_counts <- initial_df %>% dplyr::select(all_of(id_var))
  
  # Loop through the rest of the files, combine the counts by gene, and add to the dataframe
  for (i in 1:length(file_list)) {
    message("Appending ", file_list[i])
    counts <- reader(file_list[i]) %>% dplyr::select(all_of(c(reads_column, id_var)) )
    if(sorted){combined_counts <- cbind(combined_counts, counts[, reads_column]) }
    if(!sorted){
      warning("Sorted files makes this much, much faster. This is not an error but encouragement.")
      combined_counts <- merge(combined_counts, counts, by=id_var)}
    colnames(combined_counts)[ncol(combined_counts)] <- basename(file_list[i])
    }
  
  # Keep only distinct ids
  combined_counts <- combined_counts %>% distinct(get(id_var), .keep_all = T) %>% dplyr::select(-`get(id_var)`)
  # Replace any NA values with 0
  combined_counts[is.na(combined_counts)] <- 0
  
  # Return the combined dataframe
  return(combined_counts)
  }

change_df_datatype <- function(df, fxn) {
  require(dplyr)
  rows <- rownames(df)
  df <- df %>% mutate_all(fxn)
  rownames(df) <- rows
  return(df)
}

###########################################################################
#
#                                 CODE
#
###########################################################################
# Add code here
#
# Read in metadata and counts
metadata <- read_xls('data/GSE65705/GSE65705_RNAseq_ClinicalVars_GEO_13Nov2014.xls')
counts <- collect_counts_from_folder("data/GSE65705/GSE65705_RAW", reads_column="tag_count", id_var="accession")

# Change accession number to gene IDs
counts$accession <- str_extract(counts$accession,"^[A-Z]{2}_\\d+")
counts$symbols <- mapIds(org.Hs.eg.db, keys = counts$accession, column = "SYMBOL", keytype = "ACCNUM")

# Just get one read per gene.and remove smaples that don't have metadata
counts <- counts %>% distinct(symbols, .keep_all = T) %>% 
  na.omit %>% remove_rownames %>% column_to_rownames("symbols") %>%
  dplyr::select(-accession, -contains("NonMI")) %>%
  change_df_datatype(as.numeric)

# set rownames for metadata and set rownames for counts
metadata <- metadata %>% column_to_rownames("Study ID")
colnames(metadata) <- gsub(" ", "_", colnames(metadata))
colnames(counts) <- rownames(metadata)

# make SE
gse65705_se <- make_se(counts, metadata)

# Make and run DESeq
gse65705_dds <- DESeqDataSet(gse65705_se, ~ Batch + MI_Type)
gse65705_dds <- DESeq(gse65705_dds)





###########################################################################
#
#                       Export for python
#
###########################################################################
# Get press
press <- read.csv('data/clean/press_genes.csv', header = F) %>% unlist


# Normalize counts
countsNormalized <- counts(gse65705_dds, normalize=T) %>% 
  as.data.frame %>% 
  add_missing_rows(press)
countsNormalized <- countsNormalized[press,]

# setting labels as {0:NSTEMI, 1:STEMI}
labels <- metadata$MI_Type %>% mapvalues(c("NSTEMI", "STEMI"), c(0, 1))

# We are missing like 50 press genes
write.csv(countsNormalized, quote = F, file = paste0(outdir, "GSE65705_press_counts.csv") )
write.csv(labels, quote = F, file = paste0(outdir, "GSE65705_labels.csv"), row.names = F)




