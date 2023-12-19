#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Negative Control Parameters
params.OUTDIR = "${launchDir}/results"
params.BATCH_SIZE = 400
params.VERBOSITY = 1
params.HDF5_OUTPUT = "results.hdf5"
params.PERMUTATION_HDF5_OUTPUT = "permutation_results.hdf5"
params.PERMUTATION_JSON_OUTPUT = "NO_JSON_OUTPUT"
params.MAX_PERMUTATION_TESTS = ""
params.PERMUTATION_ORDERS = "1"
params.RNG = 123
params.MAF_MATCHING_RELTOL = 0.05
params.N_RANDOM_VARIANTS = 10
params.TMLE_SAVE_EVERY = 100
params.ARROW_OUTPUT = "dataset.arrow"
params.ESTIMATOR_KEY = "TMLE"

include { longest_prefix } from './utils.nf'
include { EstimationWorkflow } from './estimation.nf'

process GeneratePermutationTestsData {
    container "olivierlabayle/negative-controls:cvtmle"
    publishDir "${params.OUTDIR}/permutation_tests", mode: 'symlink', pattern: '*.arrow'
    publishDir "${params.OUTDIR}/permutation_tests/estimands", mode: 'symlink', pattern: '*.jls'
    label "bigmem"
    
    input:
        path dataset
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.jls", emit: estimands

    script:
        limit = params.MAX_PERMUTATION_TESTS == "" ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_permutation_data.jl \
        ${dataset} ${results} \
        ${limit} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --orders=${params.PERMUTATION_ORDERS} \
        --chunksize=${params.BATCH_SIZE} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

process GenerateRandomVariantsTestsData {
    container "olivierlabayle/negative-controls:cvtmle"
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    
    input:
        path trans_actors
        path bgenfiles
        path results

    output:
        path "random_variants_estimands.jls"

    script:
        trans_actors_prefix = longest_prefix(trans_actors)
        bgen_prefix = longest_prefix(bgenfiles)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_random_variant_parameters.jl \
        ${trans_actors_prefix} ${results} ${bgen_prefix} \
        --out=random_variants_estimands.jls \
        --p=${params.N_RANDOM_VARIANTS} \
        --reltol=${params.MAF_MATCHING_RELTOL} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

workflow {
    results_file = Channel.value(file("${params.OUTDIR}/${params.HDF5_OUTPUT}"))

    // Permutation Tests
    dataset = Channel.value(file("${params.OUTDIR}/${params.ARROW_OUTPUT}"))
    estimator_config = Channel.value(file("${params.ESTIMATOR_FILE}"))
    keep_ic = false
    do_svp  = false
    pval_threshold = params.PVAL_THRESHOLD
    save_every = params.TMLE_SAVE_EVERY
    hdf5_output = params.PERMUTATION_HDF5_OUTPUT
    json_output = params.PERMUTATION_JSON_OUTPUT

    GeneratePermutationTestsData(
        dataset, 
        results_file
    )
    EstimationWorkflow(
        GeneratePermutationTestsData.output.dataset,
        GeneratePermutationTestsData.output.estimands.flatten(),
        estimator_config,
        keep_ic,
        do_svp,
        pval_threshold,
        save_every,
        hdf5_output,
        json_output
    )

    // Random Variants parameter files generation
    if (params.STUDY_DESIGN == "FROM_ACTORS") {
        bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS", checkIfExists: true).collect()
        GenerateRandomVariantsTestsData(trans_actors, bgen_files, results_file)
    }
}