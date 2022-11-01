process SieveVarianceEstimation {
    container "olivierlabayle/tl-core:0.2.0"
    publishDir "$params.OUTDIR/hdf5files", mode: 'symlink', pattern: "*.hdf5"
    publishDir "$params.OUTDIR/summaries", mode: 'symlink', pattern: "*.csv"

    input:
        path tmle_file
        path GRM_ids
        path GRM_matrix

    output:
        path "${outprefix}.csv", emit: summary
        path "${outprefix}.hdf5", emit: hdf5file
    
    script:
        outprefix = tmle_file.getName().replace("hdf5", "sieve")
        """julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/sieve_variance.jl \
        $tmle_file GRM $outprefix \
        --nb-estimators=$params.NB_VAR_ESTIMATORS --max-tau=$params.MAX_TAU --pval=$params.PVAL_SIEVE
        """
}