include { TMLE; MergeOutputs } from '../modules/estimation.nf'
include { GenerateSummaryPlots } from '../modules/plot.nf'

workflow EstimationWorkflow {
    take:
        dataset
        estimands_configs
        estimators_config
        keep_ic
        do_svp
        pval_threshold
        save_every
        hdf5_output
        json_output

    main:
        // Run the estimation process for each estimands configuration
        TMLE(
            dataset,
            estimands_configs,
            estimators_config,
            keep_ic,
            do_svp,
            pval_threshold,
            save_every
        )
        // Merge results files together
        MergeOutputs(TMLE.out.collect(), hdf5_output, json_output)

        // Generate Plots
        GenerateSummaryPlots(MergeOutputs.out.hdf5_file)

    emit:
        hdf5_result = MergeOutputs.out.hdf5_file
}