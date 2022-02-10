process VariantRun {
    container "olivierlabayle/tmle-epistasis:0.3.0"
    label "bigmem"
    publishDir "$params.OUTDIR/estimates/", mode: 'symlink'

    input:
        path bgenfiles
        path phenotypefile
        path confoundersfile
        path estimatorfile
        path queryfile
        path phenotypelist_file
    
    output:
        path "*.{hdf5, jls}"
    
    script:
        def phenotype_list = phenotypelist_file.getName() != 'NONE' ? "--phenotypes-list $phenotypelist_file" : ''
        def adaptive_cv = params.ADAPTIVE_CV == true ? '--adaptive-cv' : ''
        def save_full = params.SAVE_FULL == true ? '--save-full' : ''
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/ukbb.jl $phenotypefile $confoundersfile $queryfile $estimatorfile $phenotype_list $adaptive_cv $save_full"
    
}
