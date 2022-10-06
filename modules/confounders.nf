process FlashPCA {
    cpus 8
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
    container "olivierlabayle/tl-core:new_strategies"
    publishDir "$params.OUTDIR/covariates/", mode: 'symlink'
    memory '16 GB'
    
    input:
        path flashpca_out
    
    output:
        path "pcs.csv"
    
    script:
        "julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/prepare_confounders.jl --input $flashpca_out --output pcs.csv adapt"
}
