#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.PHENOTYPES_LIST = "NONE"


process filterASB {
    container "docker://olivierlabayle/ukbb-estimation-pipeline:0.1.0"

    input:
        path asb_snp_files

    output: 
        path "filtered_asb_snps.csv", emit: filtered_asb_snps

    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/filter_asb.jl --out filtered_asb_snps.csv ${asb_snp_files.join(" ")}"
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
        julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/generate_queries.jl $filtered_asb_snps $trans_actors -o queries -s $bgen_sample -t $params.THRESHOLD $exclude
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
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/prepare_phenotypes.jl $binary_phenotypes $continuous_phenotypes $bridge phenotypes.csv --withdrawal-list $withdrawal_list"
}


process TMLE {
    container "olivierlabayle/tmle-epistasis:0.1.4"
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


include { generateCovariates } from './modules/covariates.nf'


workflow {
    // Generate queries
    asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
    trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
    excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
    bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)

    filterASB(asb_snp_ch.collect())
    generateQueries(filterASB.out.filtered_asb_snps, trans_actors, bgen_files_ch.toList(), excluded_snps)

    // generate covariates
    qc_file = Channel.value(file("$params.QC_FILE"))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    generateCovariates(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)

    // generate phenotypes
    binary_phenotypes = Channel.value(file("$params.BINARY_PHENOTYPES"))
    continuous_phenotypes = Channel.value(file("$params.CONTINUOUS_PHENOTYPES"))
    bridge = Channel.value(file("$params.GENEATLAS_BRIDGE"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))


    generatePhenotypes(binary_phenotypes, continuous_phenotypes, bridge, withdrawal_list)

    // generate TMLE arguments tuples
    // Binary phenotypes:
    binary_phen = Channel.value(file("$params.BINARY_PHENOTYPES", checkIfExists: true))
                        .splitCsv(sep: "\s", limit: 1)
                        .flatten()
                        .filter { it != "FID" & it != "IID"}
    if (params.PHENOTYPES_LIST != "NONE") { 
        binary_phen = binary_phen.filter(params.PHENOTYPES_LIST)
    }
    binary_est_phen = binary_phen.combine(Channel.fromPath("$params.BINARY_ESTIMATORFILE", checkIfExists: true))

    // Binary phenotypes:
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
    TMLE(bgen_files_ch.collect(), generatePhenotypes.out, generateCovariates.out, phenotypes_estimators_queries)
}