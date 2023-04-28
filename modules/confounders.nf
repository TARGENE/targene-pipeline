process FlashPCA {
    label "multithreaded"
    container "ktetleycampbell/flashpca:1.0"

    input:
        path bedfiles
    
    output:
        path "pcs.txt"
    
    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile $prefix --ndim $params.NB_PCS --numthreads $task.cpus"
}

process AdaptFlashPCA {
    container "olivierlabayle/tl-core:up_tmle_dep"
    publishDir "$params.OUTDIR/covariates/", mode: 'symlink'
    label 'bigmem'
    
    input:
        path flashpca_out
    
    output:
        path "pcs.csv"
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input $flashpca_out --output pcs.csv adapt
        """
}
