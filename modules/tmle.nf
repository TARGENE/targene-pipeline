process TMLE {
    container "olivierlabayle/tmle-epistasis:0.1.5"
    label "bigmem"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path bgenfiles
        path phenotypefile
        path confoundersfile
        path estimatorfile
        tuple file(phenotypelist_file), file(queryfile)
    
    output:
        path "estimates*"
    
    script:
        outfile = "estimates_" + queryfile.baseName + "_" + phenotype + ".csv"
        "julia --project=/GenesInteraction.jl --startup-file=no /GenesInteraction.jl/ukbb_epistasis.jl $phenotypefile $confoundersfile $queryfile $estimatorfile $outfile --phenotypes-list $phenotypelist_file"
    
}