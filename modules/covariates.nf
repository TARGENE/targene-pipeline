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
    container "olivierlabayle/tl-core:v0.1.0"
    publishDir "$params.OUTDIR/covariates/", mode: 'symlink'
    memory '16 GB'
    
    input:
        path flashpca_out
    
    output:
        path "pcs.csv"
    
    script:
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/prepare_confounders.jl --input $flashpca_out --output pcs.csv adapt"
}

process MergeExtraCovariatesAndPCs{
    container "olivierlabayle/ukbmain:v0.1.0"
    publishDir "$params.OUTDIR/covariates/", mode: 'symlink'
    memory '16 GB'
    
    input:
        path csv1
        path csv2

    output:
        "covariates.csv"

    script:
        """
        julia --project=/UKBMain/ /UKBMain/bin/csvmerge.jl --startup-file=no \
        $csv1 $csv2 covariates.csv
        """
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
