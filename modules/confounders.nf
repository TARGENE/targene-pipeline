include { longest_prefix } from './utils'

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
        qc_file_flag = qcfile.getName() != 'NO_QC_FILE' ? "--qc-file ${qcfile}" : '' 
        ld_blocks_flag = ld_blocks.getName() != 'NO_LD_BLOCKS' ? "--ld-blocks-file ${ld_blocks}" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/targenecore.jl \
            filter-chromosome ${input_prefix} filtered.${input_prefix} ${traits} \
            ${qc_file_flag} \
            --maf-threshold=${params.MAF_THRESHOLD} \
            ${ld_blocks_flag}
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
        tuple val(genotypes_id), path(bed_files)
    
    output:
        tuple val(genotypes_id), path("${output_prefix}*")

    script:
        output_prefix = "${genotypes_id}.merged"
        input_prefix = longest_prefix(bed_files)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/targenecore.jl \
        merge-beds ${input_prefix} ${output_prefix}
        """

}

process SampleQCFilter {
    label 'bigmem'
    label 'plink_image'
    publishDir "${params.OUTDIR}/iid_genotypes", mode: 'symlink'

    input:
        tuple val(genotypes_id), path(merged_bed_files)
    
    output:
        tuple val(genotypes_id), path("${output_prefix}*")

    script:
        input_prefix = "${genotypes_id}.merged"
        output_prefix = "${input_prefix}.qc"
        "plink2 --bfile ${input_prefix} --make-bed --hwe 1e-10 --geno --mind --out ${output_prefix}"
}

process FlashPCA {
    label "multithreaded"
    label 'pca_image'
    publishDir "${params.OUTDIR}/pve", mode: 'symlink', pattern: "*pve*"

    input:
        tuple val(genotypes_id), path(bedfiles)
    
    output:
        tuple val(genotypes_id), path("pcs.${genotypes_id}.txt"), emit: pcs
        tuple val(genotypes_id), path("pve.${genotypes_id}.txt")
    
    script:
        input_prefix = bedfiles[0].toString().minus('.bed')
        "/home/flashpca-user/flashpca/flashpca --bfile ${input_prefix} --ndim ${params.NB_PCS} --numthreads ${task.cpus} --suffix .${genotypes_id}.txt"
}
