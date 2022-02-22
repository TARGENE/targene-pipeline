process SieveVarianceEstimation {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    label "bigmem"

    input:
        path estimate_file
        path GRM_ids
        path GRM_matrix
        val nb_estimators

    output:
        path "${estimate_file.getName()}"
    
    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/sieve_variance.jl ${estimate_file.getName()} GRM --nb-estimators=$nb_estimators"
}