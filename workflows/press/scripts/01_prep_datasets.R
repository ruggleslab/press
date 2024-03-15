###########################################################################
#
#                            prep_datasets
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: prep_datasets

#" The JSON object contains three key-value pairs:
#" * "normalization": This key corresponds to the method used for normalization. In this case, "mor" is used.
#" * "geneset": This key corresponds to the path of the selected features file. In this case, the path is "output/feature_selection/rfe/logreg/selected_features.csv".
#" * "outdir": This key corresponds to the output directory where the results will be stored. In this case, the output directory is "output/reduce_press/rfeLOGREG_mor/".
# Load in the json from argparse
parser <- argparse::ArgumentParser()
parser$add_argument("--json", help = "path to the config file")
args <- parser$parse_args()
params <- jsonlite::fromJSON(args$json)

# Output directory:
experiment <- "datasets"
outdir <- file.path(params$outdir, paste0(experiment))
dir.create(outdir, showWarnings = F)

# Would be nice to have one script that is going to test variious normalization and transformation methods

#======================== LIBRARIES ========================#
library(tidyverse)
library(readxl)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
suppressMessages(source("https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R"))

#" Prepares the data for learning.
#"
#" This function takes in count data, metadata, an output directory, and a set of parameters.
#" It normalizes the count data, selects a gene set, and prepares the data for DESeq2.
#"
#" @param counts A data frame of count data.
#" @param meta A data frame of metadata.
#" @param outdir A string specifying the output directory.
#" @param params A list of parameters. This should be a list that includes the following elements:
#"   * normalization: A string specifying the normalization method. For example, "mor".
#"   * geneset: A string specifying the path to the selected features file. For example, "output/feature_selection/rfe/logreg/selected_features.csv".
#"   * outdir: A string specifying the output directory. For example, "output/reduce_press/rfeLOGREG_mor/".
#" @return list A list containing the normalized counts and the metadata.
prep_data <- function(
    counts, 
    meta, 
    outdir,
    params = params
    ) {
    # function body
    dir.create(outdir, showWarnings = F, recursive = T)
    # get the intersection of the counts and the meta
    common <- intersect(colnames(counts), rownames(meta))
    if (length(common) == 0) {
        stop("No common samples between counts and meta")
    }
    counts <- counts[, common]
    meta <- meta[common, ]

    # get the params
    normalize <- params$normalization
    geneset <- pull(read.csv(params$geneset))

    # add missing rows to the counts
    counts <- add_missing_rows(counts, geneset)

    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts[, common],
        colData = meta[common, ],
        design = ~ 1
    )
    dds <- DESeq2::estimateSizeFactors(dds)
    counts.norm <- normalize_counts(
        dds, method = normalize,
        log = ifelse(normalize %in% c("cpm", "tmm", "mor"), TRUE, FALSE)
        )
    label <- meta

    write.csv(counts.norm[, common], file.path(outdir, "normalized_counts.csv"))
    write.csv(counts[, common], file.path(outdir, "raw_counts.csv"))
    write.csv(label[common, ], file.path(outdir, "label.csv"))
    return(list(counts = counts.norm, label = label))
}

#======================== PACE ========================#
# Load data
pace_counts <- read.table("output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/concat_count_files.txt", row.names = 1)
hypercohort_metatable <- read.csv("output/hyper_geneset_creation/run18_hyper60_hypo40_AGRCONTROL/hypercohort_metatable.csv", row.names = 1)

hypercohort_metatable$hypercohort_inrnaseq_AP <- plyr::mapvalues(hypercohort_metatable$hypercohort_inrnaseq_AP, from = c("hyper", "nothyper"), to = c(1, 0))

pace_data <- prep_data(
    pace_counts, 
    hypercohort_metatable %>% drop_na(hypercohort_inrnaseq_AP), 
    file.path(outdir, "derivation"),
    params = params
    )

#======================== DUKE ========================#
duke_counts <- read.csv("data/duke_validation_run3/dukerawcounttable_conv.csv", row.names = 1)
duke_subjects_id <- rownames(read.csv("data/clean/duke_longitudinal_group.csv", row.names = 1))
duke_metadata <- read.csv("data/duke_validation_run3/dukemetatable_sel.csv", row.names = 1)
duke_metadata <- duke_metadata %>% 
    filter(characteristic__subject_id %in% duke_subjects_id) %>%
    mutate(
        hypercohort = case_when(
            characteristic__epi_max_05 > 60 ~ 1,
            characteristic__epi_max_05 < 40 ~ 0,
            TRUE ~ NA_real_
        )
    )
# with(duke_metadata, table(hypercohort, cohort))

duke_data <- prep_data(
    duke_counts, 
    duke_metadata %>% filter(cohort == "group1") %>% select(hypercohort),
    file.path(outdir, "duke_t1"),
    params = params
    )
duke_data <- prep_data(
    duke_counts, 
    duke_metadata %>% filter(cohort == "group2") %>% select(hypercohort),
    file.path(outdir, "duke_t2"),
    params = params
    )
duke_all <- prep_data(
    duke_counts, 
    duke_metadata %>% select(hypercohort),
    file.path(outdir, "duke_all"),
    params = params
    )

#======================== HARP ========================#
harp_rawcounts <- read.csv("data/harp_plt/harp_allmicontrol__rawcounttable.csv", header = TRUE, row.names = 1)
harp_metadata <- read.csv("data/harp_plt/harp_allmicontrol__metatable.csv", header = TRUE, row.names = 1)
harp_metadata_addon <- read.csv("data/harp_plt/harp-platelet-metadata-ADDON-20220823.csv", header = TRUE, row.names = 1)
harp_metadata_addon2 <- read_xlsx("data/harp_plt/Metadata HARP for WB RNAseq 3_14_2022x.xlsx")

# subset the metadata addon to the harp metadata
harp_metadata_addon <- harp_metadata_addon[rownames(harp_metadata),]
harp_metadata <- harp_metadata_addon
rownames(harp_metadata) <- gsub("-", "\\.", rownames(harp_metadata))
harp_metadata_addon2$Subject.ID <- gsub("-", "\\.", harp_metadata_addon2$Subject.ID)

# list from Tessa Barrett of who to include
harp_selection <- c(
# MI-CAD   
"HARP-01-0246-1",
"HARP-01-0098-1",
"HARP-01-0052-1",
"HARP-01-0039-1",
"HARP-01-0033-1",
"HARP-01-0035-1",
"HARP-01-109-1",
"HARP-01-0204-1",
"HARP-01-0034-1",
"HARP-01-0135-1",
"HARP-01-0177-1",
"HARP-01-0245-1",
"HARP-01-0236-1",
"HARP-01-0070-1",
"HARP-01-0043-1",
"HARP-01-0055-1",
"HARP-01-0057-1",
"HARP-01-0174-1",
"HARP-01-0241-1",
"HARP-01-0043-1",
# Obstructive
"HARP-01-5061-1",
"HARP-01-5062-1",
"HARP-01-5057-1",
"HARP-01-5012-2",
"HARP-01-5013-1",
"HARP-01-5024-1",
"HARP-01-5027-1",
"HARP-01-5033-1",
"HARP-01-5046-1"
) %>% 
gsub("-", "\\.", .)

# get a list of the IDs we are missing from the metadata
harp_selection[!harp_selection %in% rownames(harp_metadata)]

# change HARP.01.0135.1 to MI-CAD
harp_metadata[rownames(harp_metadata) == "HARP.01.0135.1", "Angiography.Report"] <- "MI-CAD"

# merge the four groups into just MI and Control
harp_metadata <- harp_metadata %>%
    filter(rownames(.) %in% harp_selection) %>%
    mutate(
        MI_v_Ctrl = ifelse(Angiography.Report == "MI-CAD", "MI-CAD",
        ifelse(Angiography.Report != "MINOCA", "Control", NA)),
        MI_v_Obstr = ifelse(Angiography.Report == "MI-CAD", "MI-CAD",
        ifelse(Angiography.Report == "Obstructive", "Obstructive", NA))
    )
harp_metadata$MI_v_Ctrl <- plyr::mapvalues(harp_metadata$MI_v_Ctrl, from = c("MI-CAD", "Control"), to = c(1, 0))
harp_metadata$MI_v_Obstr <- plyr::mapvalues(harp_metadata$MI_v_Obstr, from = c("MI-CAD", "Obstructive"), to = c(1, 0))

# add the bmi
harp_metadata$bmi <- harp_metadata_addon2 %>%
  as.data.frame() %>%
  column_to_rownames("Subject.ID") %>%
  .[rownames(harp_metadata), "BMI"]

harp_data <- prep_data(
    harp_rawcounts,
    harp_metadata %>% select(MI_v_Ctrl),
    file.path(outdir, "harp"),
    params = params
    )

#======================== SLE ========================#
sle_counts <- read.table("data/raw_validation/platelet-sle/rna_processing/filtrawcounttab.txt", header = TRUE, row.names = 1)
sle_meta <- read.table("./data/raw_validation/platelet-sle/rna_processing/metatable_filt.txt", sep = "\t", header = TRUE, row.names = 1)
rownames(sle_meta) <- gsub("-", "\\.", paste0("X", rownames(sle_meta)))

sle_meta$Diagnosis <- plyr::mapvalues(sle_meta$Diagnosis, from = c("sle", "control"), to = c(1, 0))

sle_data <- prep_data(
    sle_counts, 
    sle_meta, 
    file.path(outdir, "sle"),
    params = params
    )

#======================== COVID ========================#
covid_counts <- read.table("data/raw_validation/covid-platelet/rna_processing/filtrawcounttab.txt", header = TRUE, row.names = 1)
covid_meta <- read.table("data/raw_validation/covid-platelet/rna_processing/metatable_filt.txt", sep = "\t", header = TRUE, row.names = 1)
rownames(covid_meta) <- gsub("-", "\\.", rownames(covid_meta))

covid_meta$Cohort <- plyr::mapvalues(covid_meta$Cohort, from = c("COVID", "Control"), to = c(1, 0))

covid_data <- prep_data(
    covid_counts, 
    covid_meta %>% select(Cohort), 
    file.path(outdir, "covid"),
    params = params
)

#======================== PAD ========================#
pace_dds <- readRDS("data/pace/dds.rds")
pace_counts <- assay(pace_dds)
pace_metatable <- as.data.frame(colData(pace_dds))
pace_data <- prep_data(
    pace_counts, 
    pace_metatable %>% select(PAD),
    file.path(outdir, "PAD"),
    params = params
)

#======================== MACLE2 ========================#
pace_dds <- readRDS("data/pace/dds.rds")
pace_counts <- assay(pace_dds)
pace_metatable <- as.data.frame(colData(pace_dds))
# remove the NA values
pace_metatable <- pace_metatable %>% drop_na(censor_MACLE2)
pace_counts <- pace_counts[, rownames(pace_metatable)]
pace_data <- prep_data(
    pace_counts,
    pace_metatable %>% select(censor_MACLE2),
    file.path(outdir, "MACLE2"),
    params = params
)

#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
