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

