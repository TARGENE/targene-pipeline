process VariantRun {
    container "olivierlabayle/tmle-epistasis:0.22.0"
    label "bigmem"

    input:
        val mode
        path bgenfiles
        path phenotypefile
        path confoundersfile
        path estimatorfile
        tuple file(queryfile), file(phenotypelist_file)
    
    output:
        path "${mode}_estimates.csv"
    
    script:
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/ukbb.jl $phenotypefile $confoundersfile $queryfile $estimatorfile ${mode}_estimates.csv --phenotypes-list $phenotypelist_file $mode"
    
}
