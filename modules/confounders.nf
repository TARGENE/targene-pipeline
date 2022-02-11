process buildPCs {
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

process adapt_flashpca {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    publishDir "$params.OUTDIR/confounders", mode: 'symlink'
    memory '16 GB'
    
    input:
        path flashpca_out
    
    output:
        path "confounders.csv"
    
    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/prepare_confounders.jl --input $flashpca_out --output confounders.csv adapt"
}


workflow generatePCs{
    take:
        iid_genotypes
    main:
        pcs = buildPCs(iid_genotypes)
        adapt_flashpca(pcs)
    emit:
        adapt_flashpca.out
}