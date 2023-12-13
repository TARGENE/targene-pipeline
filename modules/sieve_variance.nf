process SieveVarianceEstimation {
    container "olivierlabayle/targeted-estimation:cv_tmle"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path tmle_files
        path GRM_ids
        path GRM_matrix

    output:
        path "svp.hdf5"
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no /opt/bin/tmle sieve-variance-plateau \
        tmle_result \
        --n-estimators=$params.NB_SVP_ESTIMATORS \
        --max-tau=$params.MAX_SVP_THRESHOLD \
        --estimator-key=$params.SVP_ESTIMATOR_KEY \
        --verbosity=$params.VERBOSITY
        """
}

process MergeOutputs {
    container "olivierlabayle/targeted-estimation:cv_tmle"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path tmle_files
        val hdf5_output
        val json_output

    output:
        path "${hdf5_output}"
        path "${json_output}", optional: true

    script:
        json_output = json_output != "NO_JSON_OUTPUT" ? "--outputs.json.filename=${json_output}" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no /opt/bin/tmle make-summary \
        tmle_result \
        --outputs.hdf5.filename=${params.HDF5_OUTPUT} \
        ${json_output}
        """
}