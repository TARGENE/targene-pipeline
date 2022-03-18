#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.PHENOTYPES_LIST = "NO_FILE"
params.QUERIES_MODE = "given"
params.THRESHOLD = 0.9
params.ADAPTIVE_CV = true
params.SAVE_FULL = false

include { IIDGenotypes } from './modules/genotypes.nf'
include { generatePCs } from './modules/confounders.nf'
include { queriesFromASBxTransActors; filterASB; queriesFromQueryFiles } from './modules/queries.nf'
include { phenotypesFromGeneAtlas as BridgeContinuous; phenotypesFromGeneAtlas as BridgeBinary } from './modules/phenotypes.nf'
include { TMLE as TMLEContinuous; TMLE as TMLEBinary} from './modules/tmle.nf'
include { GRMPart; AggregateGRM } from './modules/grm.nf'
include { SieveVarianceEstimation } from './modules/sieve_variance.nf'


workflow generateIIDGenotypes {
    qc_file = Channel.value(file("$params.QC_FILE"))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)

    emit:
        IIDGenotypes.out
}

workflow generateGRM {
    take:
        iid_genotypes

    main:
        grm_parts = Channel.from( 1..params.GRM_NSPLITS )
        GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)

        // Aggregate files
        AggregateGRM(GRMPart.out.collect())

    emit:
        grm_ids = AggregateGRM.out.grm_ids
        grm_matrix = AggregateGRM.out.grm_matrix

}


workflow generateConfounders {
    take:
        iid_genotypes

    main:
        generatePCs(iid_genotypes.collect())

    emit:
        generatePCs.out
}


workflow generateQueries{
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
    binary_phenotypes = Channel.value(file("$params.BINARY_PHENOTYPES"))
    continuous_phenotypes = Channel.value(file("$params.CONTINUOUS_PHENOTYPES"))
    bridge = Channel.value(file("$params.GENEATLAS_BRIDGE"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))

    BridgeContinuous(continuous_phenotypes, bridge, withdrawal_list, params.PHENOTYPES_LIST, "continuous_phenotypes.csv")

    BridgeBinary(binary_phenotypes, bridge, withdrawal_list, params.PHENOTYPES_LIST, "binary_phenotypes.csv")
    
    emit:
        continuous = BridgeContinuous.out
        binary = BridgeBinary.out
}


workflow generateEstimates {
    take:
        continuous_phenotypes_file
        binary_phenotypes_file
        queries_files
        confounders_file

    main:
        bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))

        // compute TMLE estimates
        TMLEContinuous(bgen_files_ch.collect(), continuous_phenotypes_file, confounders_file, estimator_file, queries_files, "Real")
        TMLEBinary(bgen_files_ch.collect(), binary_phenotypes_file, confounders_file, estimator_file, queries_files, "Bool")

    emit:
        continuous = TMLEContinuous.out
        binary = TMLEBinary.out
}


workflow generateVarianceEstimates {
    take:
        estimate_file
        GRM_ids
        GRM_matrix
    
    main:
        SieveVarianceEstimation(estimate_file, GRM_ids, GRM_matrix, params.NB_VAR_ESTIMATORS)
    
    emit:
        SieveVarianceEstimation.out
}


workflow {
    // Generate queries
    generateQueries()

    // Generate IID Genotypes
    generateIIDGenotypes()

    // generate confounders
    generateConfounders(generateIIDGenotypes.out)

    // generate GRM
    generateGRM(generateIIDGenotypes.out)

    // generate phenotypes
    generatePhenotypes()

    // generate estimates
    generateEstimates(
        generatePhenotypes.continuous.out,
        generatePhenotypes.binary.out,
        generateQueries.out.flatten(), 
        generateConfounders.out
    )
    
    // generate variance estimates
    //generateVarianceEstimates(generateEstimates.out, generateGRM.out.grm_ids, generateGRM.out.grm_matrix)
}