process filterBED {
    label 'bigmem'
    label 'targenecore_image'
    publishDir "${params.OUTDIR}/qc_filtered_chromosomes", mode: 'symlink'

    input:
        tuple val(chr_id), file(bedfiles)
        path qcfile
        path ld_blocks
        path traits

    output:
        path "filtered.*", emit: filtered_bedfiles

    script:
        input_prefix = bedfiles[0].toString().minus('.bed')
        qc_file = qcfile.getName() != 'NO_QC_FILE' ? "--qc-file ${qcfile}" : '' 
        ld_blocks = ld_blocks.getName() != 'NO_LD_BLOCKS' ? "--ld-blocks-file ${ld_blocks}" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/targenecore.jl \
            filter-chromosome ${input_prefix} filtered.${input_prefix} ${traits} \
            ${qc_file} \
            --maf-threshold=${params.MAF_THRESHOLD} \
            ${ld_blocks}
        """
}


process thinByLD {
    label 'bigmem'
    label 'plink_image'
    publishDir "${params.OUTDIR}/ld_pruned_chromosomes", mode: 'symlink'

    input:
        path flashpca_excl_reg
        path bedfiles

    output:
        path "LDpruned.*"

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        """
        plink2 --memory ${task.memory.toMega()} --bfile ${prefix} --indep-pairwise 1000 50 0.05 --exclude range ${flashpca_excl_reg}
        plink2 --memory ${task.memory.toMega()} --bfile ${prefix} --extract plink2.prune.in --make-bed --out LDpruned.${prefix}
        """
}


process mergeBEDS{
    label 'bigmem'
    label 'targenecore_image'
    publishDir "${params.OUTDIR}/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        path "ukbb_merged*"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/targenecore.jl \
        merge-beds LDpruned. ukbb_merged
        """

}

process SampleQCFilter {
    label 'bigmem'
    label 'plink_image'
    publishDir "${params.OUTDIR}/iid_genotypes", mode: 'symlink'

    input:
        path merged_bed_files
    
    output:
        path "qc_filtered*"

    script:
        "plink2 --bfile ukbb_merged --make-bed --hwe 1e-10 --geno --mind --out qc_filtered"
}

process FlashPCA {
    label "multithreaded"
    label 'pca_image'

    input:
        path bedfiles
    
    output:
        path "pcs.txt"
    
    script:
        prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile ${prefix} --ndim ${params.NB_PCS} --numthreads ${task.cpus}"
}

process AdaptFlashPCA {
    publishDir "${params.OUTDIR}/covariates/", mode: 'symlink'
    label 'bigmem'
    label 'targenecore_image'

    input:
        path flashpca_out
    
    output:
        path "pcs.csv"
    
    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/targenecore.jl \
        ${flashpca_out} pcs.csv
        """
}