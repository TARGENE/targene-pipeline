include { TMLE; GenerateOutputs } from '../modules/estimation.nf'

workflow EstimationWorkflow {
    take:
        dataset_and_estimands
        estimators_config

    main:
        // Run the estimation process for each estimands configuration
        TMLE(
            dataset_and_estimands,
            estimators_config,
        )

        // Generate TarGene Outputs
        GenerateOutputs(TMLE.out.collect())

    emit:
        hdf5_result = GenerateOutputs.out
}