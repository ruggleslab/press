###########################################################################
#
#                            prep_datasets
#
###########################################################################
# Author: Matthew Muller
# Date: 2024-02-25
# Script Name: prep_datasets

# Load in the json from argparse
parser <- argparse::ArgumentParser()
parser$add_argument('--json', help = 'path to the config file')
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
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/general_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/plotting_functions.R')
source('https://raw.githubusercontent.com/mattmuller0/Rtools/main/stats_functions.R')

dataset_dir <- glue("{params$outdir}/datasets")
modeling_dir <- glue("{params$outdir}/modeling")
post_dir <- glue("{params$outdir}/postprocessing")
dir.create(post_dir, showWarnings = F)

# so we want to take the modeling results and compile them into a more readable format
# I think we should take the results and make a heatmap and a report for each model
# we can calculate the ORs and the p-values and add them as annotations there.

#======================== Modeling Overview ========================#
dir.create(file.path(post_dir, "modeling_overview"), showWarnings = F)
# get the summary file from the modeling/evalutation directory
summary_file <- list.files(modeling_dir, pattern = "summary.csv", full.names = T, recursive = T)
summary_df <- read.csv(summary_file, row.names = 1)

prediction_files <- list.files(modeling_dir, pattern = "predictions.csv", full.names = T, recursive = T)
predictions <- lapply(prediction_files, read.csv)
names(predictions) <- stringr::str_split_i(string = prediction_files, pattern = "/", n = -2)

## let's get the odds ratios and p-values for each model's predictive scores in the PACE data
pace_meta <- read.csv(glue("{dataset_dir}/PAD/label.csv"), row.names = 1)
pace_predictions <- read.csv(glue("{modeling_dir}/evaluation/PAD/predictions.csv"), row.names = 1)
paceDat <- cbind(pace_meta, pace_predictions)
write.csv(paceDat, file.path(post_dir, "pace_predictions.csv"))

#======================== Events ========================#
dir.create(file.path(post_dir, "events"), showWarnings = F)
censors <- grep("censor_", colnames(paceDat), value = T)
hr_table <- hazard_ratios_table(
  paceDat, 
  'predictions', censors, 
  per_sd = T,
  censor_prefix = 'censor_'
  time_prefix = 'time_to_'
  )
write.csv(hr_table, file.path(post_dir, "events", "hazard_ratios_per_sd.csv"))

hr_table$HR_ci_upper[hr_table$HR_ci_upper > 5] <- 5
forestPlot <- ggplot(hr_table, aes(x = estimate, y = censor, color = `p.value` < 0.05)) +
    geom_pointrange(aes(xmin = HR_ci_lower, xmax = HR_ci_upper)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    theme_matt(16) +
    lims(x = c(0, 5)) +
    theme(legend.position = 'bottom') +
    labs(x = 'Hazard Ratio Per SD', y = NULL, title = NULL, color = 'p < 0.05')
ggsave(file.path(outdir, 'events_by_press', 'hazard_ratios.pdf'), forestPlot, width = 12, height = 8)

#======================== LTA ========================#
dir.create(file.path(post_dir, "lta"), showWarnings = F)
lta <- grep('^(epi_|col_|adp_|ser_|aa_*)(?!.*slope.*|.*lag.*|.*mfi.*|.*mlty.*)', colnames(paceDat), value = T, perl = T)
lta_or_table <- purrr::map_df(
    vars, 
    ~glm(
        na.omit(stats_df[, .x]) ~ scores, 
        data = stats_df %>% drop_na(.x, labels)
        ) %>% 
        broom::tidy() %>% 
        filter(term == 'predictions') %>%
        mutate(term = .x, conf.up = estimate + 1.96 * std.error, conf.low = estimate - 1.96 * std.error) %>%
        dplyr::select(term, estimate, conf.low, conf.up, p = p.value)
    )
write.csv(lta_or_table, file.path(post_dir, "lta", "lta_or_table.csv"))
lta_or_plot <- ggplot(lta_or_table, aes(x = estimate, y = term, color = p < 0.05)) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.up)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    theme_matt(16) +
    lims(x = c(0, 5)) +
    theme(legend.position = 'bottom') +
    labs(x = 'Odds Ratio', y = NULL, title = NULL, color = 'p < 0.05')

#======================== END ========================#
writeLines(capture.output(sessionInfo()), file.path(outdir, "session.log"))
save.image(file.path(outdir, "image.RData"))
