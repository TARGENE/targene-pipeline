#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"

process filterASB {
    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"

    input:
        path asb_snp_files

    output: 
        path "filtered_asb_snps.csv", emit: filtered_asb_snps

    script:
        "julia --project --startup-file=no bin/filter_asb.jl --out filtered_asb_snps.csv ${asb_snp_files.join(" ")}"
}


process generateQueries {

    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"

    input:
        path filtered_asb_snps
        path trans_actors
        path chr_files
        path excluded_snps

    output:
        path "queries/*.toml", emit: queries

    script:
        bgen_sample = chr_files.find { it.toString().endsWith("bgen") }
        def exclude = excluded_snps.name != 'NO_FILE' ? "--exclude $excluded_snps" : ''
        """
        mkdir -p queries
        julia --project --startup-file=no bin/generate_queries.jl $filtered_asb_snps $trans_actors -o queries -s $bgen_sample -t $params.THRESHOLD $exclude
        """
}

process generatePhenotypes {
    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"
    label "bigmem"

    input:
        path binary_phenotypes
        path continuous_phenotypes
        path bridge
        path withdrawal_list
    
    output:
        path "phenotypes.csv"
    
    script:
        "julia --project --startup-file=no bin/process_phenotypes.jl $binary_phenotypes $continuous_phenotypes $bridge phenotypes.csv --withdrawal-list $withdrawal_list"
}


process TMLE {
    container "olivierlabayle/tmle-epistasis:0.1.0"

    input:
        path phenotypefile
        path confoundersfile
        path queryfile
        path estimatorfile
    
    output:
        path "estimates.csv"
    
    script:
        "julia --project --startup-file=no ukbb_epistasis.jl $phenotypefile $confoundersfile $queryfile $estimatorfile estimates.csv --phenotype $phenotype"
    
}


include { generateCovariates } from './modules/covariates.nf'


workflow {
    // Generate queries
    // asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
    // trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
    // excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
    // bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)

    // filterASB(asb_snp_ch.collect())
    // generateQueries(filterASB.out.filtered_asb_snps, trans_actors, bgen_files_ch.toList(), excluded_snps)

    // // generate covariates
    // qc_file = Channel.value(file("$params.QC_FILE"))
    // flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    // ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
    // bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    // generateCovariates(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)

    // generate phenotypes
    binary_phenotypes = Channel.value(file("$params.BINARY_PHENOTYPES"))
    continuous_phenotypes = Channel.value(file("$params.CONTINUOUS_PHENOTYPES"))
    bridge = Channel.value(file("$params.GENEATLAS_BRIDGE"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))


    generatePhenotypes(binary_phenotypes, continuous_phenotypes, bridge, withdrawal_list)

    // compute TMLE estimates
}