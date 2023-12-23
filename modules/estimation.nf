#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process MergeOutputs {
    container "olivierlabayle/targeted-estimation:argparse"
    publishDir "$params.OUTDIR", mode: 'symlink'
    label "bigmem"
    
    input:
        path tmle_files
        val hdf5_output
        val json_output

    output:
        path "${hdf5_output}", emit: hdf5_file
        path "${json_output}", optional: true, emit: json_file

    script:
        json_option = json_output != "NO_JSON_OUTPUT" ? "--json-output=${json_output}" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --sysimage=/TargetedEstimation.jl/TMLESysimage.so --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/tmle.jl merge \
        tmle_result \
        ${json_option} \
        --hdf5-output=${hdf5_output}
        """
}

process TMLE {
    container "olivierlabayle/targeted-estimation:argparse"
    publishDir "$params.OUTDIR/tmle_outputs/", mode: 'symlink', pattern: "*.hdf5"
    label "bigmem"
    label "multithreaded"

    input:
        path data
        path estimands_file
        path estimator_file
        val keep_ic
        val do_svp
        val pval_threshold
        val save_every
    
    output:
        path "${hdf5out}"
    
    script:
        basename = "tmle_result." + estimands_file.getName().take(estimands_file.getName().lastIndexOf('.'))
        hdf5out = basename + ".hdf5"
        pval_option = keep_ic == true ? ",${pval_threshold}" : ""
        sample_ids = do_svp == true ? ",true" : ""
        output_option = "--hdf5-output=${hdf5out}${pval_option}${sample_ids}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --sysimage=/TargetedEstimation.jl/TMLESysimage.so --project=/TargetedEstimation.jl --threads=${task.cpus} --startup-file=no /TargetedEstimation.jl/tmle.jl tmle \
        $data \
        --estimands=$estimands_file \
        --estimators=$estimator_file \
        $output_option \
        --chunksize=$save_every \
        """
}

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
        TMLE(
            dataset,
            estimands_configs,
            estimators_config,
            keep_ic,
            do_svp,
            pval_threshold,
            save_every
        )
        MergeOutputs(TMLE.out.collect(), hdf5_output, json_output)

    emit:
        hdf5_result = MergeOutputs.out.hdf5_file
}

// workflow {
//     dataset = Channel.value(file("${params.AGGREGATED_DATASET}"))
//     estimands_configs = Channel.value(file("${params.ESTIMAND_CONFIG_FILES}"))
//     estimators_config = Channel.value(file("${params.ESTIMATOR_FILE}", checkIfExists: false))
//     hdf5_output = "${params.HDF5_OUTPUT}"
//     json_output = "${params.JSON_OUTPUT}"
//     EstimationWorkflow(dataset, estimands_configs, estimators_config, hdf5_output, json_output)
// }