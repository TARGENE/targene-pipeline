process SieveVarianceEstimation {
    container "olivierlabayle/targeted-estimation:0.6"
    publishDir "$params.OUTDIR/hdf5files/sieve", mode: 'symlink', pattern: "*.hdf5"
    publishDir "$params.OUTDIR/csvs", mode: 'symlink', pattern: "*.csv"

    input:
        path tmle_files
        path GRM_ids
        path GRM_matrix

    output:
        path "sieve_variance.hdf5", emit: hdf5_file
        path "sieve_variance.csv", emit: csv_file
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/scripts/sieve_variance.jl \
        tmle GRM sieve_variance \
        --nb-estimators=$params.NB_VAR_ESTIMATORS --max-tau=$params.MAX_TAU
        """
}

process MergeOutputs {
    container "olivierlabayle/targeted-estimation:0.6"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path tmle_files
        path sieve_files
        val outpath

    output:
        path "${outpath}"

    script:
        tmle_prefix = "tmle"
        sieve_prefix = sieve_files.getName() == "NO_SIEVE_FILE" ? "" : "--sieve-prefix sieve_variance"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/scripts/merge_summaries.jl \
        ${tmle_prefix} ${outpath} ${sieve_prefix}
        """
}