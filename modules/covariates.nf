process filterBED{
    label 'bigmem'

    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"

    input:
        tuple val(chr_id), file(bedfiles)
        path qcfile
        path ld_blocks

    output:
        path "filtered.*", emit: filtered_bedfiles

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "julia --project=$projectDir --startup-file=no $projectDir/bin/prepare_confounders.jl --input $prefix --output filtered.$prefix --qcfile $qcfile --maf-threshold $params.MAF_THRESHOLD --ld-blocks $ld_blocks filter"

}


process thinByLD{
    label 'bigmem'
    container "docker://olivierlabayle/plink2:0.1.0"

    input:
        path flashpca_excl_reg
        path bedfiles

    output:
        path "LDpruned.*"

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        """
        plink2 --bfile $prefix --indep-pairwise 1000 50 0.05 --exclude range $flashpca_excl_reg
        plink2 --bfile $prefix --extract plink2.prune.in --make-bed --out LDpruned.$prefix
        """
}


process mergeBEDS{
    label 'bigmem'

    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"
    
    input:
        path files
    
    output:
        path "ukbb_merged*"

    script:
        "julia --project=$projectDir --startup-file=no $projectDir/bin/prepare_confounders.jl --input LDpruned. --output ukbb_merged merge"

}

process buildPCs {
    cpus 8
    container "docker://ktetleycampbell/flashpca:latest"

    input:
        path bedfiles
    
    output:
        path "pcs.txt"
    
    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile $prefix --ndim $params.NB_PCS --numthreads $task.cpus"
}


workflow generateCovariates{
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
    main:
        filtered_bedfiles = filterBED(bed_files, qc_file, ld_blocks)
        ld_pruned = thinByLD(flashpca_excl_reg, filtered_bedfiles)
        merged = mergeBEDS(ld_pruned.collect())
        buildPCs(merged)
    emit:
        buildPCs.out
}