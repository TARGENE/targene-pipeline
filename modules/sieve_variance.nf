process SieveVarianceEstimation {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    publishDir "$params.OUTDIR/hdf5files/sieve_variances", mode: 'symlink'
    label "unlimited_vmem"
    label "multithreaded"

    input:
        tuple val(rsprefix), file(estimate_file)
        path GRM_ids
        path GRM_matrix

    output:
        tuple val(rsprefix), file("${rsprefix}_sieve_variance.hdf5")
    
    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/sieve_variance.jl $rsprefix GRM ${rsprefix}_sieve_variance.hdf5 --nb-estimators=$params.NB_VAR_ESTIMATORS --max-tau=$params.MAX_TAU --pval=$params.PVAL_SIEVE"
}