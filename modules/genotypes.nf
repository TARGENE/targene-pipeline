process filterBED{
    label 'bigmem'
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/qc_filtered_chromosomes", mode: 'symlink'

    input:
        tuple val(chr_id), file(bedfiles)
        path qcfile
        path ld_blocks
        path sample_ids

    output:
        path "filtered.*", emit: filtered_bedfiles

    script:
        prefix = bedfiles[0].toString().minus('.bed')
        """
        julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/prepare_confounders.jl \
        --input $prefix --output filtered.$prefix --qcfile $qcfile --maf-threshold $params.MAF_THRESHOLD --ld-blocks $ld_blocks --sample-ids $sample_ids filter
        """

}


process thinByLD{
    label 'bigmem'
    container "olivierlabayle/plink2:0.1.0"
    publishDir "$params.OUTDIR/ld_pruned_chromosomes", mode: 'symlink'

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
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/merged_genotypes", mode: 'symlink'
    
    input:
        path files
    
    output:
        path "ukbb_merged*"

    script:
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/prepare_confounders.jl --input LDpruned. --output ukbb_merged merge"

}

process SampleQCFilter {
    label 'bigmem'
    container "olivierlabayle/plink2:0.1.0"
    publishDir "$params.OUTDIR/iid_genotypes", mode: 'symlink'

    input:
        path merged_bed_files
    
    output:
        path "qc_filtered*"

    script:
        "plink2 --bfile ukbb_merged --make-bed --hwe 1e-10 --geno --mind --out qc_filtered"
}


workflow IIDGenotypes{
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
        sample_ids

    main:
        filtered_bedfiles = filterBED(bed_files, qc_file, ld_blocks, sample_ids)
        ld_pruned = thinByLD(flashpca_excl_reg, filtered_bedfiles)
        mergeBEDS(ld_pruned.collect())
        SampleQCFilter(mergeBEDS.out.collect())

    emit:
        SampleQCFilter.out
}