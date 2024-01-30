include { longest_prefix } from '../modules/utils.nf'
include { EstimationWorkflow } from '../modules/estimation.nf'

process GeneratePermutationTestsData {
    container "olivierlabayle/tl-core:cvtmle"
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
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint=${params.POSITIVITY_CONSTRAINT} \
        --verbosity=${params.VERBOSITY} \
        --batch-size=${params.BATCH_SIZE} \
        permutation-tests ${dataset} ${results} \
        ${limit} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --orders=${params.PERMUTATION_ORDERS} \
        --rng=${params.RNG} \
        --max-attempts=${params.MAX_PERMUTATION_ATTEMPTS} \
        """
}

process GenerateRandomVariantsTestsData {
    container "olivierlabayle/tl-core:cvtmle"
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    
    input:
        path variants_to_randomize
        path bgenfiles
        path results

    output:
        path "random_variants_estimands.jls"

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_random_variant_estimands.jl \
        ${variants_to_randomize} ${results} ${bgen_prefix} \
        --out=random_variants_estimands.jls \
        --p=${params.N_RANDOM_VARIANTS} \
        --reltol=${params.MAF_MATCHING_RELTOL} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

workflow PERMUTATION_TEST {
    // Permutation Tests
    results_file = Channel.value(file("${params.RESULTS_FILE}"))
    dataset = Channel.value(file("${params.AGGREGATED_DATASET}"))
    estimator_config = Channel.value(file("${params.ESTIMATOR_FILE}"))
    keep_ic = false
    do_svp  = false
    pval_threshold = params.PVAL_THRESHOLD
    save_every = params.TMLE_SAVE_EVERY
    hdf5_output = params.HDF5_OUTPUT
    json_output = params.JSON_OUTPUT

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
}

workflow RANDOMIZATION_TEST {
    results_file = Channel.value(file("${params.RESULTS_FILE}"))
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
    variants_to_randomize = Channel.value(file("$params.VARIANTS_TO_RANDOMIZE", checkIfExists: true))
    GenerateRandomVariantsTestsData(variants_to_randomize, bgen_files, results_file)
}