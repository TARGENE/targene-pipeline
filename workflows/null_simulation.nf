include { EXTRACT_TRAITS } from './traits.nf'
include { PCA } from './pca.nf'
include { EstimationInputs } from '../subworkflows/estimation_inputs.nf'
include { NullSimulationEstimation; AggregateSimulationResults } from '../modules/simulations.nf'

workflow NULL_SIMULATION {
    // DEPRECATED: LEGACY TO BE REMOVED
    bqtls_file = Channel.value(file("$params.BQTLS"))
    transactors_files = Channel.fromPath("$params.TRANS_ACTORS").collect()
    extra_confounders = Channel.value(file("$params.EXTRA_CONFOUNDERS"))
    extra_treatments = Channel.value(file("$params.ENVIRONMENTALS"))
    extra_covariates = Channel.value(file("$params.EXTRA_COVARIATES"))

    // Workflow specific channels
    estimators = Channel.fromPath(params.ESTIMATORS, checkIfExists: true)
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    sample_sizes = Channel.fromList(params.SAMPLE_SIZES)
    rngs = Channel.fromList(params.RNGS)
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()

    // Extract Traits
    traits = EXTRACT_TRAITS()
    
    // PCA
    pcs = PCA()

    // generate main dataset and estimand configuration files
    EstimationInputs(
        bgen_files,
        traits,
        pcs,
        estimands_file,
        bqtls_file,
        transactors_files,
        extra_confounders,
        extra_treatments,
        extra_covariates,
    )

    // Null Simulation Estimation
    validated_estimands = EstimationInputs.out.estimands.flatten()
    validated_dataset = EstimationInputs.out.aggregated_dataset
    bootstrap_grid = estimators.combine(validated_estimands).combine(sample_sizes).combine(rngs)
    simulation_results = NullSimulationEstimation(validated_dataset, bootstrap_grid)

    // Aggregation of Estimation Results
    AggregateSimulationResults(simulation_results.collect(), "null_simulation_results.hdf5")
}