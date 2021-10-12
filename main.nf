#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


process filterASB {
    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"

    input:
        path asb_snp_files

    output: 
        path "filtered_asb_snps.csv", emit: filtered_asb_snps

    script:
        "julia --project=$projectDir --startup-file=no $projectDir/bin/filter_asb.jl --out filtered_asb_snps.csv ${asb_snp_files.join(" ")}"
}


process generateQueries {
    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"

    input:
        path filtered_asb_snps
        path trans_actors
        path chr_files

    output:
        path "queries/*.toml", emit: queries

    script:
        bgen_sample = chr_files.find { it.toString().endsWith("bgen") }
        """
        mkdir -p queries
        julia --project=$projectDir --startup-file=no $projectDir/bin/generate_queries.jl $filtered_asb_snps $trans_actors -o queries -s $bgen_sample -t $params.THRESHOLD
        """
}

include { generateCovariates } from './modules/covariates.nf'


workflow {
    asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
    trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)

    bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true)

    qc_file = Channel.value(file("$params.QC_FILE"))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS"))

    // Merge and filter allelic-specific binding SNPs
    filterASB(asb_snp_ch.collect())

    // Generate queries
    generateQueries(filterASB.out.filtered_asb_snps, trans_actors, bgen_files_ch.toList())

    // filterBED
    generateCovariates(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)
}