params.PERMUTATION_HDF5_OUTPUT = "permutation_results.hdf5"
params.PERMUTATION_JSON_OUTPUT = "NO_JSON_OUTPUT"

// Permutation Tests Parameters
params.MAX_PERMUTATION_TESTS = null
params.PERMUTATION_ORDERS = "1"
params.RNG = 123
params.MAF_MATCHING_RELTOL = 0.05
params.N_RANDOM_VARIANTS = 10

include { longest_prefix } from './utils.nf'

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
        limit = params.MAX_PERMUTATION_TESTS == null ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_permutation_data.jl \
        ${dataset} ${results} \
        ${limit} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
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
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

workflow NegativeControl {
    results_file = Channel.value(file("${params.OUTDIR}/${params.HDF5_OUTPUT}"))

    // Permutation Tests
    dataset = Channel.value(file("${params.OUTDIR}/${params.ARROW_OUTPUT}"))
    estimator_file = Channel.value(file("${params.ESTIMATOR_FILE}"))
    GeneratePermutationTestsData(dataset, results_file)
    TMLE(
        GeneratePermutationTestsData.output.dataset, 
        GeneratePermutationTestsData.output.estimands.flatten(), 
        estimator_file
    )
    MergeOutputs(TMLE.out.collect(), params.PERMUTATION_HDF5_OUTPUT, params.PERMUTATION_JSON_OUTPUT)
    
    // Random Variants parameter files generation
    if (params.STUDY_DESIGN == "FROM_ACTORS") {
        bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS", checkIfExists: true).collect()
        GenerateRandomVariantsTestsData(trans_actors, bgen_files, results_file)
    }
}