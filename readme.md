# PRESS Platelet Reactivity Signature Code Repository

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Description

This repository contains the code used to create PRESS. Is is mainly written in R and Python.

The code as is relies on a config file to run. The corresponding manuscript can be found [here](https://www.nature.com/articles/s41467-024-50994-7).

The PRESS Model can be run either using the 03_run_press.py script or by downloading and using the PRESS model from this repository. If you are downloading the model, we recommend you use the included scaler unless you believe there will be significant batch effects in your data. Data should be preprocessed and normalized using the log2 median of ratios method.

## Usage

* 01_press_derivation_script.R
  * This script was used to perform the analysis of the platelet hyper/hypo geneset creation and analysis
* 02_press_modeling.ipynb
  * This notebook was used to test the PRESS model
* 03_run_press.py
    * This script was used to run the PRESS model
* 04_harp_validation.R
    * This script was used to validate the PRESS model in the HARP cohort

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

* Matthew Muller
* <mm12865@nyu.edu>
* Project Link: [GitHub Repository](https://github.com/ruggleslab/press)
* Zenodo Link: [![DOI](https://zenodo.org/badge/816318277.svg)](https://zenodo.org/doi/10.5281/zenodo.12537340)
* Manuscript Link: TBD
