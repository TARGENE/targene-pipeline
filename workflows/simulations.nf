include { CreateEstimatorsConfigChannel } from '../modules/utils.nf'
include { EXTRACT_TRAITS } from './traits.nf'
include { PCA } from './pca.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { RealisticSimulationInputs; NullSimulationEstimation; RealisticSimulationEstimation; AggregateRealisticSimulationResults; AggregateNullSimulationResults; DensityEstimation } from '../modules/simulations.nf'

workflow NULL_SIMULATION {
    // Workflow specific channels
    estimators = CreateEstimatorsConfigChannel(params.ESTIMATORS_CONFIG)
    estimands_files = Channel.value(file("$params.ESTIMANDS_CONFIG"))
    sample_sizes = Channel.fromList(params.SAMPLE_SIZES)
    rngs = Channel.fromList(params.RNGS)
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect().toList()

    // PCA
    PCA()

    // generate main dataset and estimand configuration files
    pcs_and_genotypes = PCA.out.pcs.combine(bgen_files)
    EstimationInputs(
        pcs_and_genotypes,
        PCA.out.traits,
        estimands_files
    )

    // Null Simulation Estimation
    estimation_inputs = EstimationInputs.out.multiMap { dataset, estimands ->
        dataset: dataset
        estimands: estimands

    }
    validated_estimands = estimation_inputs.estimands.flatten()
    validated_dataset = estimation_inputs.dataset.collect()
    bootstrap_grid = estimators.combine(validated_estimands).combine(sample_sizes).combine(rngs)
    simulation_results = NullSimulationEstimation(validated_dataset, bootstrap_grid)

    // Aggregation of Estimation Results
    AggregateNullSimulationResults(simulation_results.collect())
}

workflow REALISTIC_SIMULATION {
    // Workflow specific channels
    ga_trait_table = Channel.value(file(params.GA_TRAIT_TABLE, checkIfExists: true))
    estimators = CreateEstimatorsConfigChannel(params.ESTIMATORS_CONFIG)
    estimands_files = Channel.value(file("$params.ESTIMANDS_CONFIG"))
    sample_sizes = Channel.fromList(params.SAMPLE_SIZES)
    rngs = Channel.fromList(params.RNGS)
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect().toList()
    
    // PCA
    PCA()

    // Realistic Simulation Inputs
    pcs_and_genotypes = PCA.out.pcs.combine(bgen_files)
    simulation_inputs = RealisticSimulationInputs(
        pcs_and_genotypes,
        estimands_files.collect(),
        PCA.out.traits,
        ga_trait_table
    )
    dataset = simulation_inputs.dataset.collect()

    // Density Estimation
    density_estimates = DensityEstimation(
        dataset,
        simulation_inputs.conditional_densities.flatten(),
    ).collect()

    // Estimation
    validated_estimands = simulation_inputs.estimands.flatten()
    bootstrap_grid = estimators.combine(validated_estimands).combine(sample_sizes).combine(rngs)
    simulation_results = RealisticSimulationEstimation(
        dataset,
        density_estimates,
        bootstrap_grid
    )

    // Aggregation of Estimation Results
    AggregateRealisticSimulationResults(simulation_results.collect(), dataset, density_estimates)
}