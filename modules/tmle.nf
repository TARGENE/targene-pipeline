process TMLE {
    container "olivierlabayle/tmle-epistasis:0.21.0"
    label "bigmem"

    input:
        path bgenfiles
        path phenotypefile
        path confoundersfile
        path estimatorfile
        tuple file(queryfile), file(phenotypelist_file)
    
    output:
        path "estimates.csv"
    
    script:
        "julia --project=/GenesInteraction.jl --startup-file=no /GenesInteraction.jl/ukbb_epistasis.jl $phenotypefile $confoundersfile $queryfile $estimatorfile estimates.csv --phenotypes-list $phenotypelist_file"
    
}