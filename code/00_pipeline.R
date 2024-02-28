library(tidyverse)

# Load in the json from argparse
parser <- argparse::ArgumentParser()
parser$add_argument('--json', help = 'path to the config file')
args <- parser$parse_args()
params <- jsonlite::fromJSON(args$json)

# make the outdir
subdir <- file.path(params$outdir, paste0('norm_', params$normalization, '__', 'gs_', stringr::str_split_i(params$geneset, '/', -2)))
dir.create(subdir, showWarnings = F, recursive = T)