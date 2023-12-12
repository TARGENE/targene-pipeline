process SieveVarianceEstimation {
    container "olivierlabayle/targeted-estimation:cv_tmle"
    publishDir "$params.OUTDIR/hdf5files/sieve", mode: 'symlink', pattern: "*.hdf5"
    publishDir "$params.OUTDIR/csvs", mode: 'symlink', pattern: "*.csv"

    input:
        path tmle_files
        path GRM_ids
        path GRM_matrix

    output:
        path "svp.hdf5"
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no tmle sieve-variance-plateau \
        result \
        --nb-estimators=$params.NB_SVP_ESTIMATORS \
        --max-tau=$params.MAX_SVP_THRESHOLD \
        --estimator_key=$params.SVP_ESTIMATOR_KEY \
        --verbosity=$params.VERBOSITY
        """
}

process MergeOutputs {
    container "olivierlabayle/targeted-estimation:cv_tmle"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path tmle_files
        val outpath

    output:
        path "${outpath}"

    script:
        json_output = params.JSON_OUTPUT != "NO_JSON_OUTPUT" ? "--outputs.json.filename=$params.JSON_OUTPUT" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no tmle make-summary \
        result \
        --outputs.hdf5.filename=${params.HDF5_OUTPUT} \
        ${json_output}
        """
}