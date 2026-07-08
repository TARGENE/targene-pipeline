include { TMLE } from '../modules/estimation.nf'

workflow EstimationWorkflow {
    take:
        dataset_and_estimands
        estimators_config

    main:
        // Run the estimation process for each estimands configuration
        tmle_results = TMLE(dataset_and_estimands.combine(estimators_config)).collect()

    emit:
        hdf5_result = tmle_results
}