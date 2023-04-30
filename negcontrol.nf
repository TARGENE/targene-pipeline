#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Global Parameters
params.OUTDIR = "${launchDir}/results"
params.RESULTS_FILE = ""
params.TMLE_INPUT_DATASET = ""
params.VERBOSITY = 1

// Permutation Tests Parameters
params.MAX_PERMUTATION_TESTS = null
params.PVAL_COL = "TMLE_PVALUE"
params.PVAL_THRESHOLD = 0.05
params.PERMUTATION_CHUNKSIZE = 1000
params.PERMUTATION_ORDERS = "1"
params.PERMUTATION_RNG = 123
params.ESTIMATORFILE = ""


include { TMLE } from './modules/tmle.nf'


process GeneratePermutationTestsData {
    container "olivierlabayle/negative-controls:initial_pipeline"
    publishDir "${params.OUTDIR}/neg_control", mode: 'symlink'
    label "bigmem"
    
    input:
        path dataset
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.bin", emit: parameters

    script:
        limit = params.MAX_PERMUTATION_TESTS == null ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_permutation_data.jl \
        ${dataset} ${results} \
        ${limit} \
        --pval-col=${params.PVAL_COL} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --orders=${params.PERMUTATION_ORDERS} \
        --chunksize=${params.PERMUTATION_CHUNKSIZE} \
        --rng=${params.PERMUTATION_RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

workflow {
    dataset = Channel.value(file("${params.TMLE_INPUT_DATASET}"))
    results_file = Channel.value(file("${params.RESULTS_FILE}"))
    estimatorfile = Channel.value(file("${params.ESTIMATORFILE}"))

    GeneratePermutationTestsData(dataset, results_file)
    TMLE(
        GeneratePermutationTestsData.output.dataset, 
        GeneratePermutationTestsData.output.parameters, 
        estimatorfile
    )
}