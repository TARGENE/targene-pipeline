process TMLE {
    container "olivierlabayle/tmle-epistasis:0.3.0"
    label "bigmem"
    label "multithreaded"

    input:
        path bgenfiles
        path phenotypefile
        path confoundersfile
        path estimatorfile
        path queryfile
        val target_type
    
    output:
        path "*.hdf5"
    
    script:
        def adaptive_cv = params.ADAPTIVE_CV == true ? '--adaptive-cv' : ''
        def save_full = params.SAVE_FULL == true ? '--save-full' : ''
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/ukbb.jl $phenotypefile $confoundersfile $queryfile $estimatorfile --target-type $target_type $adaptive_cv $save_full"
    
}
