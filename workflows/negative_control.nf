include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { GeneratePermutationTestsData; GenerateRandomVariantsTestsData} from '../modules/negative_control.nf'

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