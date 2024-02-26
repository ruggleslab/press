# PRESS Platelet Reactivity Signature Code Repository

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

This repository contains the code used to create PRESS. Is is mainly written in R and Python.

The code as is relies on a config file to run. The config file is located in the config folder and is called `params.json`.

## Table of Contents

* PRESS
 ├──code
 │  ├──data-wrangling.R
 │  ├──extract_weights.py
 │  ├──feature_selection.py
 │  ├──gene_qc_investigation.R
 │  ├──GSE65705_data_cleaning.R
 │  ├──harp_predictions.R
 │  ├──PACE_GEO_SUBM.R
 │  ├──pace_press_scoring_tiles.R
 │  ├──pace_scoring_investigation.R
 │  ├──platelet-hyper-hypo-analysis-3-agr-control.R
 │  ├──platelet-pace_analysis_master_2.R
 │  ├──run_press.py
 │  ├──scratch.html
 │  └──scratch.Rmd
 ├──config
 │  └──params.json
 ├──models
 │  ├──jobs
 │  │  ├──press451.joblib
 │  │  └──press451_with_scaler.joblib
 │  ├──press_reduction_scaler.joblib
 │  └──weights
 │     ├──press451
 │     │  └──estimators.npy
 │     └──press451_with_scaler
 │        ├──steps.npy
 │        └──voting__estimators.npy
 ├──notebooks
 │  ├──data-cleaning.ipynb
 │  ├──feature_selection.ipynb
 │  ├──harp_hyper_v_hypo_predictions.ipynb
 │  ├──model-exploring.ipynb
 │  ├──model_testing.ipynb
 │  └──model_testing_new_genesets.ipynb
 ├──readme.md
 └──src
    ├──data_functions.py
    ├──feature_selection.py
    ├──tosh_rnaseq_scripts
    │  ├──deseq_functions.R
    │  ├──geneset_analysis_functions.R
    │  ├──mgc_plotting_functions.R
    │  ├──overlap_finder_function.R
    │  ├──rnaseq_analysis_README.docx
    │  ├──rnaseq_processing_functions.R
    │  ├──ssGSEA_custom.R
    │  ├──symbol_species_conversion_functions.R
    │  └──TEMPLATE_analysis_master.R
    └──validation_functions.py

## Usage

* platelet-hyper-hypo-analysis-3-agr-control.R
  * This script was used to perform the analysis of the platelet hyper/hypo geneset creation and analysis
* model_testing.ipynb
  * This notebook was used to test the PRESS model
* run_press.py
    * This script was used to run the PRESS model

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

* Matthew Muller
* <mm12865@nyu.edu>
* Project Link: [GitHub Repository](https://github.com/mattmuller0/press)
* Manuscript Link: TBD
