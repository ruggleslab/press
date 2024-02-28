###########################################################################
#
#                                 HEADER
#
###########################################################################
# Author: Matthew Muller
# 
# Date: 2023-04-11
# 
# Script Name: PACE GEO SUBMISSION
# 
# Notes:
# This script is to manipulate data to format for the pace data submission work.

  # Output directory
experiment <- "PACE_GEO_SUBM"
outdir <- file.path("output", paste0(experiment, "/"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

###########################################################################
#
#                                 LIBRARIES
#
###########################################################################
packages <- c(
  "devtools",
  "languageserver",
  "tidyverse",
  "ggplot2",
  "devtools",
  "BiocManager",
  "SummarizedExperiment",
  "DESeq2",
  "edgeR",
  "dplyr"
)

for (pkg in packages) {library(pkg, character.only = TRUE, quietly = TRUE)}


# LOAD FUNCTIONS
# space reserved for sourcing in functions
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source_url('https://raw.githubusercontent.com/mattmuller0/Rtools/main/deseq_helper_functions.R')


###########################################################################
#
#                                 CODE
#
###########################################################################
# Load in input count tables and metadata

# Load in count tables
load('data/summarized_experiments/se_pace_pltRNA.rdata')
rawCounts <- se_pace_plt %>% assay()
press_counts <- read.csv('data/hyper_feature_outtable.csv', row.names = 1)

# Get sample names and gene names
sampleNames <- colnames(rawCounts)
geneNames <- rownames(rawCounts)

# Load in metadata
metadata <- read.csv('data/metadata_9_20220422.csv', row.names = 1)
metadata <- metadata[sampleNames, ]

hypercohort_metadata <- read.csv('data/hypercohort_metatable.csv', row.names = 1)
hypercohort_metadata <- hypercohort_metadata[sampleNames, ]

print(colnames(metadata))

metadata <- metadata %>% filter(Cohort == "PACE", timepoint == "baseline") %>% dplyr::select(
  race, age, age_cat, sex1, ethnicity, antiplatelet_therapy, 
  HTN, diabetes, smoking_status, CAD_hx, PAD_hx,
  censor_MACE, censor_MALE, censor_MACLE2,
  time_to_MACE, time_to_MALE, time_to_MACLE2
  ) %>% mutate(
    epi_04um_300s_n = hypercohort_metadata$epi_04um_300s_n,
    hypercohort_binary = hypercohort_metadata$hypercohort,
    sex = sex1, sex1 = NULL,
    censor_MACLE = censor_MACLE2, censor_MACLE2 = NULL,
    time_to_MACLE = time_to_MACLE2, time_to_MACLE = NULL
    )

# Write out metadata
write.csv(metadata, file.path(outdir, 'metadata.csv'), row.names = T, quote = F)

# Write out count tables
write.csv(rawCounts, file.path(outdir, 'rawCounts.csv'), row.names = T, quote = F)
write.csv(press_counts, file.path(outdir, 'press_counts.csv'), row.names = T, quote = F)

# Print out dim of count tables and metadata
print(dim(rawCounts))
print(dim(press_counts))
print(dim(metadata))
 

# Create a list of sample names to use to grab the correct
# R1 and R2 fastq files for each sample

# bigpurple fastq address:
# /gpfs/data/ruggleslab/genomics/berger/pace_platelet/pace_platelet_v1_20200201/fastq
# bigpurple fastq fmt example:
# PACE_117_S81_L002_R1_001.fastq.gz

# Get sample names
sampleNames <- colnames(rawCounts)

# S regex is 2 digits (S[0-9]{2})
# L regex is 3 digits (L00[1-4])
# R1 and R2 regex is 1 digit
# Get R1 and R2 fastq file names
path_soft_link = '/gpfs/data/ruggleslab/genomics/berger/pace_platelet/pace_platelet_v1_20200201/fastq/'
path_real = '/gpfs/data/sequence/results/bergerlab/2019-09-27/fastq/'
pattern = paste0(sampleNames, '_S*_L*_R*.fastq.gz')
# Set which path to use
path = path_soft_link

sample_path_names <- list.files(path = path, pattern = pattern, full.names = T)
regex_sample_names <- paste0(path, pattern)

# Write out sample path names
write.table(regex_sample_names, file.path(outdir, 'regex_sample_names.txt'), row.names = F, quote = F, col.names = F)

