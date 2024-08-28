process MergeOutputs {
    publishDir "$params.OUTDIR", mode: 'symlink'
    label "bigmem"
    label 'tmle_image'

    input:
        path tmle_files

    output:
        path "${params.HDF5_OUTPUT}", emit: hdf5_file
        path "${params.JSON_OUTPUT}", optional: true, emit: json_file

    script:
        json_option = params.JSON_OUTPUT != "NO_JSON_OUTPUT" ? "--json-output=${params.JSON_OUTPUT}" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --sysimage=/TMLECLI.jl/TMLESysimage.so --project=/TMLECLI.jl --startup-file=no /TMLECLI.jl/tmle.jl merge \
        tmle_result \
        ${json_option} \
        --hdf5-output=${params.HDF5_OUTPUT}
        """
}

process TMLE {
    publishDir "$params.OUTDIR/tmle_outputs/", mode: 'symlink', pattern: "*.hdf5"
    label 'tmle_image'

    input:
        tuple path(dataset), path(estimands_file)
        path estimator_file
    
    output:
        path "${hdf5out}"
    
    script:
        basename = "tmle_result." + estimands_file.getName().take(estimands_file.getName().lastIndexOf('.'))
        hdf5out = basename + ".hdf5"
        pvalue_threhsold = params.KEEP_IC == true ? "--pvalue-threshold=${params.PVAL_THRESHOLD}" : ""
        save_sample_ids = params.SVP == true ? "--save-sample-ids" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --sysimage=/TMLECLI.jl/TMLESysimage.so --project=/TMLECLI.jl --threads=${task.cpus} --startup-file=no /TMLECLI.jl/tmle.jl tmle \
        ${dataset} \
        --estimands=${estimands_file} \
        --estimators=${estimator_file} \
        --hdf5-output=${hdf5out} \
        ${pvalue_threhsold} \
        ${save_sample_ids} \
        --chunksize=${params.TL_SAVE_EVERY}
        """
}

