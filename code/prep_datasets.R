###########################################################################
#
#                            prep_datasets
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: prep_datasets
# Output directory:
experiment <- "prep_datasets"
outdir <- file.path('output', paste0(experiment))
dir.create(outdir, showWarnings = F)

# Would be nice to have one script that is going to test variious normalization and transformation methods

#======================== LIBRARIES ========================#
library(tidyverse)

# LOAD FUNCTIONS
# space reserved for sourcing in functions
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')

params <- jsonlite::fromJSON('config/params.json')

# function to prep the data for learning
prep_data <- function(
    counts, 
    meta, 
    outdir,
    params = params
    ) {
    dir.create(outdir, showWarnings = F, recursive = T)
    # get the intersection of the counts and the meta
    common <- intersect(colnames(counts), rownames(meta))
    if (length(common) == 0) {
        stop('No common samples between counts and meta')
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
    counts.norm <- normalize_counts(dds, method = normalize)
    label <- meta

    write.csv(counts.norm[, common], file.path(outdir, 'normalized_counts.csv'))
    write.csv(counts[, common], file.path(outdir, 'raw_counts.csv'))
    write.csv(label[common, ], file.path(outdir, 'label.csv'))
    return(list(counts = counts.norm, label = label))
}

#======================== PACE ========================#
# Load data
pace_counts <- read.table('output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/concat_count_files.txt', row.names = 1)
hypercohort_metatable <- read.csv('output/hyper_geneset_creation/run18_hyper60_hypo40_AGRCONTROL/hypercohort_metatable.csv', row.names = 1)
pace_data <- prep_data(
    pace_counts, 
    hypercohort_metatable %>% drop_na(hypercohort_inrnaseq_AP), 
    file.path(outdir, 'derivation'),
    params = params
    )

pace_counts <- read.table('output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/concat_count_files.txt', row.names = 1)
hypercohort_metatable <- read.csv('output/run7_rmoutliers2_agesexcontrol_withpaired_20220422/rna_processing/metatable_in.csv', row.names = 1)
pace_data <- prep_data(
    pace_counts, 
    hypercohort_metatable, 
    file.path(outdir, 'pace'),
    params = params
    )

#======================== DUKE ========================#
duke_counts <- read.csv('data/duke_validation_run3/dukerawcounttable_conv.csv', row.names = 1)
duke_subjects_id <- rownames(read.csv('data/clean/duke_longitudinal_group.csv', row.names = 1))
duke_metadata <- read.csv('data/duke_validation_run3/dukemetatable_sel.csv', row.names = 1)
duke_metadata <- duke_metadata %>% 
    filter(characteristic__subject_id %in% duke_subjects) %>%
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
    duke_metadata %>% filter(cohort == 'group1'),
    file.path(outdir, 'duke_t1'),
    params = params
    )
duke_data <- prep_data(
    duke_counts, 
    duke_metadata %>% filter(cohort == 'group2'),
    file.path(outdir, 'duke_t2'),
    params = params
    )
duke_all <- prep_data(
    duke_counts, 
    duke_metadata,
    file.path(outdir, 'duke_all'),
    params = params
    )

#======================== HARP ========================#
harp_rawcounts <- read.csv('data/harp_plt/harp_allmicontrol__rawcounttable.csv', header = TRUE, row.names = 1)
harp_metadata <- read.csv('data/harp_plt/harp_allmicontrol__metatable.csv', header = TRUE, row.names = 1)
harp_metadata_addon <- read.csv('data/harp_plt/harp-platelet-metadata-ADDON-20220823.csv', header = TRUE, row.names = 1)
harp_metadata_addon2 <- read_xlsx('data/harp_plt/Metadata HARP for WB RNAseq 3_14_2022x.xlsx')

# subset the metadata addon to the harp metadata
harp_metadata_addon <- harp_metadata_addon[rownames(harp_metadata),]
harp_metadata <- harp_metadata_addon
rownames(harp_metadata) <- gsub("-", "\\.", rownames(harp_metadata))
harp_metadata_addon2$Subject.ID <- gsub("-", "\\.", harp_metadata_addon2$Subject.ID)

# list from Tessa Barrett of who to include
harp_selection <- c(
# MI-CAD   
'HARP-01-0246-1',
'HARP-01-0098-1',
'HARP-01-0052-1',
'HARP-01-0039-1',
'HARP-01-0033-1',
'HARP-01-0035-1',
'HARP-01-109-1',
'HARP-01-0204-1',
'HARP-01-0034-1',
'HARP-01-0135-1',
'HARP-01-0177-1',
'HARP-01-0245-1',
'HARP-01-0236-1',
'HARP-01-0070-1',
'HARP-01-0043-1',
'HARP-01-0055-1',
'HARP-01-0057-1',
'HARP-01-0174-1',
'HARP-01-0241-1',
'HARP-01-0043-1',
# Obstructive
'HARP-01-5061-1',
'HARP-01-5062-1',
'HARP-01-5057-1',
'HARP-01-5012-2',
'HARP-01-5013-1',
'HARP-01-5024-1',
'HARP-01-5027-1',
'HARP-01-5033-1',
'HARP-01-5046-1'
) %>% 
gsub("-", "\\.", .)

# get a list of the IDs we are missing from the metadata
harp_selection[!harp_selection %in% rownames(harp_metadata)]

# change HARP.01.0135.1 to MI-CAD
harp_metadata[rownames(harp_metadata) == 'HARP.01.0135.1', 'Angiography.Report'] <- 'MI-CAD'

# merge the four groups into just MI and Control
harp_metadata <- harp_metadata %>%
    filter(rownames(.) %in% harp_selection) %>%
    mutate(
        MI_v_Ctrl = ifelse(Angiography.Report == 'MI-CAD', 'MI-CAD',
        ifelse(Angiography.Report != 'MINOCA', 'Control', NA)),
        MI_v_Obstr = ifelse(Angiography.Report == 'MI-CAD', 'MI-CAD',
        ifelse(Angiography.Report == 'Obstructive', 'Obstructive', NA))
    )
harp_metadata$MI_v_Ctrl <- factor(harp_metadata$MI_v_Ctrl, levels = c('Control', 'MI-CAD'))
harp_metadata$MI_v_Obstr <- factor(harp_metadata$MI_v_Obstr, levels = c('Obstructive', 'MI-CAD'))

# add the bmi
harp_metadata$bmi <- harp_metadata_addon2 %>%
  as.data.frame() %>%
  column_to_rownames('Subject.ID') %>%
  .[rownames(harp_metadata), 'BMI']

harp_data <- prep_data(
    harp_rawcounts, 
    harp_metadata, 
    file.path(outdir, 'harp'),
    params = params
    )

#======================== SLE ========================#




#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
