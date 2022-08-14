process SieveVarianceEstimation {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/hdf5files/sieve_variances", mode: 'symlink'

    input:
        tuple val(treatment_id), file(tmle_files)
        path GRM_ids
        path GRM_matrix

    output:
        tuple val(treatment_id), file("final.${treatment_id}.sieve_variance.hdf5")
    
    script:
        """julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/sieve_variance.jl \
        final GRM final.${treatment_id}.sieve_variance.hdf5 \
        --nb-estimators=$params.NB_VAR_ESTIMATORS --max-tau=$params.MAX_TAU --pval=$params.PVAL_SIEVE
        """
}

process NoSieveEstimation {
    input:
        tuple val(treatment_id), file(tmle_files)

    output:
        tuple val(treatment_id), val("NOFILE")
        
    script:
        "exit 0"
}