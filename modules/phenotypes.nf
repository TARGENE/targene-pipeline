process phenotypesFromGeneAtlas {
    container "olivierlabayle/ukbb-estimation-pipeline:0.2.0"
    label "bigmem"

    input:
        path binary_phenotypes
        path continuous_phenotypes
        path bridge
        path withdrawal_list
    
    output:
        path "phenotypes.csv"
    
    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/prepare_phenotypes.jl $binary_phenotypes $continuous_phenotypes $bridge phenotypes.csv --withdrawal-list $withdrawal_list"
}