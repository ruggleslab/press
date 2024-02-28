#!/usr/bin/env nextflow run
// Set logging
NXF_LOG_FILE="/logs/nextflow.log"

// Define some parameters
params.normalization = 'tmm'
params.geneset_path = 'output/feature_selection/rfe/rf/selected_features.csv'
params.outdir = '/TEST/'

// Define Nextflow process for running Python script
process runCode {
    publishDir "${params.outdir}", mode: 'copy'
    script:
    """
    mkdir -p ${params.outdir}

    Rscript steps/nf_prep_datasets.R \
        --normalization ${params.normalization} \
        --geneset_path ${params.geneset_path} \
        --outdir ${params.outdir}

    python steps/nf_press_model.py \
        --geneset_path ${params.geneset_path} \
        --outdir ${params.outdir}

    echo "normalization=${params.normalization}" > ${params.outdir}/params.txt
    echo "geneset_path=${params.geneset_path}" >> ${params.outdir}/params.txt
    echo "outdir=${params.outdir}" >> ${params.outdir}/params.txt
    """
}


// Define workflow dependencies
workflow {
    // Define the order of execution
    runCode()
}