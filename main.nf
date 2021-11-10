#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.PHENOTYPES_LIST = "NONE"
params.QUERIES_MODE = "given"
params.THRESHOLD = 0.9


workflow generateConfounders {
    include { generatePCs } from './modules/confounders.nf'

    qc_file = Channel.value(file("$params.QC_FILE"))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    
    generatePCs(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)

    emit:
        generatePCs.out
}



workflow generateQueries{
    include { queriesFromASBxTransActors; filterASB; queriesFromQueryFiles } from './modules/queries.nf'

    if (params.QUERIES_MODE == "ASBxTransActors") {
        asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
        excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
        bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
        filterASB(asb_snp_ch.collect())
        queries = queriesFromASBxTransActors(filterASB.out.filtered_asb_snps, 
                                             trans_actors, 
                                             bgen_files_ch.toList(), 
                                             excluded_snps)
    }
    else if (params.QUERIES_MODE == "given"){
        queries = queriesFromQueryFiles(Channel.fromPath("$params.QUERY_FILES", checkIfExists: true))
    }
    emit:
        queries

}

workflow generatePhenotypes {
    include { phenotypesFromGeneAtlas } from './modules/phenotypes.nf'

    binary_phenotypes = Channel.value(file("$params.BINARY_PHENOTYPES"))
    continuous_phenotypes = Channel.value(file("$params.CONTINUOUS_PHENOTYPES"))
    bridge = Channel.value(file("$params.GENEATLAS_BRIDGE"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))

    phenotypesFromGeneAtlas(binary_phenotypes, continuous_phenotypes, bridge, withdrawal_list)
    
    emit:
        phenotypesFromGeneAtlas.out
}


process TMLE {
    container "olivierlabayle/tmle-epistasis:0.1.5"
    label "bigmem"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path bgenfiles
        path phenotypefile
        path confoundersfile
        tuple val(phenotype), file(estimatorfile), file(queryfile)
    
    output:
        path "estimates*"
    
    script:
        outfile = "estimates_" + queryfile.baseName + "_" + phenotype + ".csv"
        "julia --project=/GenesInteraction.jl --startup-file=no /GenesInteraction.jl/ukbb_epistasis.jl $phenotypefile $confoundersfile $queryfile $estimatorfile $outfile --phenotype $phenotype"
    
}


workflow {
    // Generate queries
    generateQueries()
    
    // generate confounders
    generateConfounders()

    // generate phenotypes
    generatePhenotypes()

    // generate TMLE arguments tuples
    bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
    // Binary phenotypes:
    binary_phen = Channel.value(file("$params.BINARY_PHENOTYPES", checkIfExists: true))
                        .splitCsv(sep: "\s", limit: 1)
                        .flatten()
                        .filter { it != "FID" & it != "IID"}
    if (params.PHENOTYPES_LIST != "NONE") { 
        binary_phen = binary_phen.filter(params.PHENOTYPES_LIST)
    }
    binary_est_phen = binary_phen.combine(Channel.fromPath("$params.BINARY_ESTIMATORFILE", checkIfExists: true))

    // Continuous phenotypes:
    continuous_phen = Channel.value(file("$params.CONTINUOUS_PHENOTYPES"))
                        .splitCsv(sep: "\s", limit: 1)
                        .flatten()
                        .filter { it != "FID" & it != "IID"}
    if (params.PHENOTYPES_LIST != "NONE") { 
        continuous_phen = continuous_phen.filter(params.PHENOTYPES_LIST)
    }
    continuous_est_phen = continuous_phen.combine(Channel.fromPath("$params.CONTINUOUS_ESTIMATORFILE", checkIfExists: true))
    // Full argument queue
    phenotypes_estimators_queries = binary_est_phen.concat(continuous_est_phen)
                                        .combine(generateQueries.out.flatten())

    // compute TMLE estimates
    TMLE(bgen_files_ch.collect(), generatePhenotypes.out, generateConfounders.out, phenotypes_estimators_queries)
}