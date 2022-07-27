#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.DECRYPTED_DATASET = "NO_FILE"
params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.QUERIES_MODE = "given"
params.CALL_THRESHOLD = 0.9
params.MINOR_CAT_FREQUENCY = 0.001
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
include { FlashPCA; AdaptFlashPCA } from './modules/confounders.nf'
include { FromASBxTransActors; FromGivenQueries } from './modules/queries.nf'
include { UKBFieldsList; UKBConv; TraitsFromUKB } from './modules/ukb_traits.nf'
include { TMLE as TMLEContinuous; TMLE as TMLEBinary} from './modules/tmle.nf'
include { PhenotypesBatches as ContinuousPhenotypesBatches; PhenotypesBatches as BinaryPhenotypesBatches; FinalizeTMLEInputs} from './modules/tmle.nf'
include { GRMPart; AggregateGRM } from './modules/grm.nf'
include { SieveVarianceEstimation } from './modules/sieve_variance.nf'
include { Summary } from './modules/summary.nf'

workflow extractTraits {
    traits_config = Channel.value(file("$params.TRAITS_CONFIG"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))
    if (params.DECRYPTED_DATASET == "NO_FILE") {
        encrypted_dataset = Channel.value(file("$params.ENCRYPTED_DATASET"))
        encoding_file = Channel.value(file("$params.ENCODING_FILE"))
        UKBFieldsList(traits_config)
        decrypted_dataset = UKBConv(UKBFieldsList.out, encrypted_dataset, encoding_file)
    }
    else {
        decrypted_dataset = Channel.value(file("$params.DECRYPTED_DATASET"))
    }

    TraitsFromUKB(decrypted_dataset, traits_config, withdrawal_list)

    emit:
        sample_ids = TraitsFromUKB.out.sample_ids
        binary_phenotypes = TraitsFromUKB.out.binary_phenotypes
        continuous_phenotypes = TraitsFromUKB.out.continuous_phenotypes
        confounders = TraitsFromUKB.out.confounders
        treatments = TraitsFromUKB.out.treatments
        covariates = TraitsFromUKB.out.covariates
}

workflow generateIIDGenotypes {
    take:
        sample_ids

    main:
        qc_file = Channel.value(file("$params.QC_FILE"))
        flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
        ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
        bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

        IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file, sample_ids)

    emit:
        IIDGenotypes.out
}

workflow geneticConfounders {
    take:
        iid_genotypes

    main:
        FlashPCA(iid_genotypes)
        AdaptFlashPCA(FlashPCA.out)

    emit:
        AdaptFlashPCA.out

}

workflow generateQueriesAndGenotypes{
    take:
        sample_ids

    main:
        bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
        excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
        if (params.QUERIES_MODE == "ASBxTransActors") {
            asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
            trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
            outputs = FromASBxTransActors(bgen_files_ch.collect(),
                                                asb_snp_ch.collect(), 
                                                trans_actors, 
                                                excluded_snps,
                                                sample_ids)
        }
        else if (params.QUERIES_MODE == "given"){
            query_files = Channel.fromPath("$params.QUERY_FILES", checkIfExists: true).collect()
            outputs = FromGivenQueries(bgen_files_ch.collect(), query_files, excluded_snps, sample_ids)
        }

    emit:
        genotypes = outputs.genotypes
        queries = outputs.queries

}

workflow generateEstimates {
    take:
        genotypes
        queries
        continuous_phenotypes
        binary_phenotypes
        genetic_confounders
        extra_confounders
        covariates
        extra_treatments

    main:
        tmle_data = FinalizeTMLEInputs(
            continuous_phenotypes,
            binary_phenotypes,
            genotypes,
            genetic_confounders,
            extra_confounders,
            covariates,
            extra_treatments
        )
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))

        // compute TMLE estimates for continuous targets
        ContinuousPhenotypesBatches(tmle_data.continuous_phenotypes)
        queries_to_continuous_phenotype_batches = queries.combine(ContinuousPhenotypesBatches.out.flatten())
        //TMLEContinuous(tmle_data.genotypes, tmle_data.continuous_phenotypes, covariates_file, estimator_file, queries_to_continuous_phenotype_batches, "Real")
        
        // compute TMLE estimates for binary targets
        BinaryPhenotypesBatches(tmle_data.binary_phenotypes)
        queries_to_binary_phenotype_batches = queries.combine(BinaryPhenotypesBatches.out.flatten())
        //TMLEBinary(genotypes_file, binary_phenotypes_file, covariates_file, estimator_file, queries_to_binary_phenotype_batches, "Bool")

        // hdf5_files = TMLEContinuous.out.flatten()
        //                 .concat(TMLEBinary.out.flatten())
        //                 .map { it -> [it.getName().split("_batch")[0], it]}
        //                 .groupTuple()

    // emit:
    //     hdf5_files
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
        Summary(tmle_files.join(sieve_files))
}

workflow {
    // Extract traits
    extractTraits()

    // Generate queries
    generateQueriesAndGenotypes(extractTraits.out.sample_ids)

    // Generate IID Genotypes
    generateIIDGenotypes(extractTraits.out.sample_ids)

    // Genetic confounders
    geneticConfounders(generateIIDGenotypes.out)

    // generate estimates
    generateEstimates(
        generateQueriesAndGenotypes.out.genotypes,
        generateQueriesAndGenotypes.out.queries.flatten(),
        extractTraits.out.continuous_phenotypes,
        extractTraits.out.binary_phenotypes,
        geneticConfounders.out,
        extractTraits.out.confounders.ifEmpty("NO_FILE"),
        extractTraits.out.covariates.ifEmpty("NO_FILE"),
        extractTraits.out.treatments.ifEmpty("NO_FILE")
    )

    // generate sieve estimates
    //generateSieveEstimates(generateEstimates.out, generateIIDGenotypes.out)

    // generate Summaries
    //generateSummaries(generateEstimates.out, generateSieveEstimates.out)
}