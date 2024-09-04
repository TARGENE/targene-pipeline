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

process GenerateOutputs {
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    label 'targenecore_image'

    input:
        path results_file

    output:
        path "*.png"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
        make-outputs tmle_result \
        --output-prefix="." \
        --verbosity=${params.VERBOSITY}
        """
}