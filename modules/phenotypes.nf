process phenotypesFromGeneAtlas {
    publishDir "$params.OUTDIR/phenotypes", mode: 'symlink'
    container "olivierlabayle/tl-core:v0.1.0"
    label "bigmem"

    input:
        path phenotype_file
        path bridge
        path withdrawal_list
        path phenotypes_list
        val outname
    
    output:
        path outname
    
    script:
        def phen_list = phenotypes_list.getName() != 'NO_FILE' ? "--phenotypes-list $phenotypes_list" : ''
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/prepare_phenotypes.jl $phenotype_file $bridge $outname --withdrawal-list $withdrawal_list $phen_list"
}