#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.PHENOTYPES_LIST = "NO_FILE"
params.QUERIES_MODE = "given"
params.THRESHOLD = 0.9
params.ADAPTIVE_CV = true
params.SAVE_FULL = false
params.PHENOTYPES_BATCH_SIZE = 1
params.GRM_NSPLITS = 100
params.MAF_THRESHOLD = 0.01
params.NB_PCS = 6
params.NB_VAR_ESTIMATORS = 0
params.MAX_TAU = 0.8
params.PVAL_SIEVE = 0.05
params.OUTDIR = "$launchDir/results"

include { IIDGenotypes } from './modules/genotypes.nf'
include { generatePCs } from './modules/confounders.nf'
include { queriesFromASBxTransActors; filterASB; queriesFromQueryFiles } from './modules/queries.nf'
include { phenotypesFromGeneAtlas as BridgeContinuous; phenotypesFromGeneAtlas as BridgeBinary } from './modules/phenotypes.nf'
include { TMLE as TMLEContinuous; TMLE as TMLEBinary} from './modules/tmle.nf'
include {PhenotypesBatches as ContinuousPhenotypesBatches; PhenotypesBatches as BinaryPhenotypesBatches} from './modules/tmle.nf'
include { GRMPart; AggregateGRM } from './modules/grm.nf'
include { SieveVarianceEstimation } from './modules/sieve_variance.nf'
include { Summary } from './modules/summary.nf'


workflow generateIIDGenotypes {
    qc_file = Channel.value(file("$params.QC_FILE"))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)

    emit:
        IIDGenotypes.out
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
    binary_phenotypes = Channel.fromPath("$params.BINARY_PHENOTYPES")
    continuous_phenotypes = Channel.fromPath("$params.CONTINUOUS_PHENOTYPES")
    bridge = Channel.fromPath("$params.GENEATLAS_BRIDGE")
    withdrawal_list = Channel.fromPath("$params.WITHDRAWAL_LIST")
    phenotypes_list = Channel.fromPath("$params.PHENOTYPES_LIST")

    BridgeContinuous(continuous_phenotypes, bridge, withdrawal_list, phenotypes_list, "continuous_phenotypes.csv")

    BridgeBinary(binary_phenotypes, bridge, withdrawal_list, phenotypes_list, "binary_phenotypes.csv")
    
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

        // compute TMLE estimates for continuous targets
        ContinuousPhenotypesBatches(continuous_phenotypes_file)
        queries_to_continuous_phenotype_batches = queries_files.combine(ContinuousPhenotypesBatches.out.flatten())
        TMLEContinuous(bgen_files_ch.collect(), continuous_phenotypes_file, confounders_file, estimator_file, queries_to_continuous_phenotype_batches, "Real")
        
        // compute TMLE estimates for binary targets
        BinaryPhenotypesBatches(binary_phenotypes_file)
        queries_to_binary_phenotype_batches = queries_files.combine(BinaryPhenotypesBatches.out.flatten())
        TMLEBinary(bgen_files_ch.collect(), binary_phenotypes_file, confounders_file, estimator_file, queries_to_binary_phenotype_batches, "Bool")

        hdf5_files = TMLEContinuous.out.flatten()
                        .concat(TMLEBinary.out.flatten())
                        .map { it -> [it.getName().split("_batch")[0], it]}
                        .groupTuple()

    emit:
        hdf5_files
}


workflow generateSieveEstimates {
    take:
        snps_tmle_files
        iid_genotypes
    
    main:
        if (params.NB_VAR_ESTIMATORS != 0){
            // Build the GRM
            grm_parts = Channel.from( 1..params.GRM_NSPLITS )
            GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)
            AggregateGRM(GRMPart.out.collect())
            // Sieve estimation
            sieve_estimates = SieveVarianceEstimation(snps_tmle_files, AggregateGRM.out.grm_ids, AggregateGRM.out.grm_matrix)
        }
        else {
            sieve_estimates = snps_tmle_files.map(it -> [it[0], "NO_FILE"])
        }
    emit:
        sieve_estimates
}

workflow generateSummaries {
    take:
        tmle_files
        sieve_files
    
    main:
        // joining on the prefix which corresponds to a tuple of SNPS
        all = tmle_files.join(sieve_files)
        Summary(all)
}

workflow {
    // Generate queries
    generateQueries()

    // Generate IID Genotypes
    generateIIDGenotypes()

    // generate confounders
    generateConfounders(generateIIDGenotypes.out)

    // generate phenotypes
    generatePhenotypes()

    // generate estimates
    generateEstimates(
        generatePhenotypes.out.continuous.first(),
        generatePhenotypes.out.binary.first(),
        generateQueries.out.flatten(), 
        generateConfounders.out.first()
    )

    // generate sieve estimates
    generateSieveEstimates(generateEstimates.out, generateIIDGenotypes.out)

    // generate Summaries
    generateSummaries(generateEstimates.out, generateSieveEstimates.out)
}