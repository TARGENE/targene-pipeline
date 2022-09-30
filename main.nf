#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.DECRYPTED_DATASET = "NO_FILE"
params.MODE = "GivenParameters"
params.PARAMETER_FILES = "NO_PARAMETER_FILE"
params.CALL_THRESHOLD = 0.9
params.POSITIVITY_CONSTRAINT = 0.01
params.SAVE_MODELS = false
params.SAVE_IC = true
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
include { UKBFieldsList; UKBConv; TraitsFromUKB } from './modules/ukb_traits.nf'
include { TMLE as TMLEContinuous; TMLE as TMLEBinary} from './modules/tmle.nf'
include { TMLEInputsFromGivenParams; TMLEInputsFromASBTrans } from './modules/tmle.nf'
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
        traits = TraitsFromUKB.out
}

workflow generateIIDGenotypes {
    take:
        traits

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

workflow generateEstimates {
    take:
        continuous_phenotypes
        binary_phenotypes
        genetic_confounders
        extra_confounders
        covariates
        extra_treatments

    main:
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))
        bgen_files = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true).collect()
        parameter_files = Channel.fromPath("$params.PARAMETER_FILES").collect()

        if (params.MODE == "ASBxTransActors") {
            asbs = Channel.fromPath("$params.ASB_FILES", checkIfExists: true).collect()
            trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
            tmle_inputs = TMLEInputsFromASBTrans(
                bgen_files,
                binary_phenotypes,
                continuous_phenotypes,
                genetic_confounders,
                extra_confounders,
                extra_treatments,
                covariates,
                asbs, 
                trans_actors, 
                parameter_files)
        }
        else if (params.MODE == "GivenParameters"){
            tmle_inputs = TMLEInputsFromGivenParams(
                bgen_files,
                binary_phenotypes,
                continuous_phenotypes,
                genetic_confounders,
                extra_confounders,
                extra_treatments,
                covariates,
                parameter_files)
        }
        // compute TMLE estimates for continuous targets
        TMLEContinuous(
            tmle_inputs.treatments,
            tmle_inputs.continuous_phenotypes, 
            tmle_inputs.confounders,
            tmle_inputs.continuous_parameters.flatten(),
            estimator_file,
            tmle_inputs.covariates.ifEmpty(file("NO_COVARIATE")),
            "Real")
        
        // compute TMLE estimates for binary targets
        TMLEBinary(
            tmle_inputs.treatments,
            tmle_inputs.binary_phenotypes, 
            tmle_inputs.confounders,
            tmle_inputs.binary_parameters.flatten(),
            estimator_file,
            tmle_inputs.covariates.ifEmpty(file("NO_COVARIATE")),
            "Bool")

        hdf5_files = TMLEContinuous.out.flatten()
                        .concat(TMLEBinary.out.flatten())
                        .map { it -> [it.getName().tokenize(".")[2], it]}
                        .groupTuple()

    emit:
        hdf5_files
}


workflow generateSieveEstimates {
    take:
        tmle_files
        iid_genotypes
    
    main:
        if (params.NB_VAR_ESTIMATORS != 0){
            // Build the GRM
            grm_parts = Channel.from( 1..params.GRM_NSPLITS )
            GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)
            AggregateGRM(GRMPart.out.collect())
            // Sieve estimation
            sieve_estimates = SieveVarianceEstimation(tmle_files, AggregateGRM.out.grm_ids, AggregateGRM.out.grm_matrix)
        }
        else {
            sieve_estimates = tmle_files.map(it -> [it[0], "NO_FILE"])
        }
    emit:
        sieve_estimates
}

workflow generateSummaries {
    take:
        tmle_files
        sieve_files
    
    main:
        // joining on the prefix which corresponds to a tuple of Treatments
        Summary(tmle_files.join(sieve_files))
}

workflow {
    // Extract traits
    extractTraits()

    // Generate IID Genotypes
    generateIIDGenotypes(extractTraits.out)

    // Genetic confounders
    geneticConfounders(generateIIDGenotypes.out)

    // generate estimates
    generateEstimates(
        extractTraits.out.continuous_phenotypes,
        extractTraits.out.binary_phenotypes,
        geneticConfounders.out,
        extractTraits.out.confounders.ifEmpty(file("NO_EXTRA_CONFOUNDER")),
        extractTraits.out.covariates.ifEmpty(file("NO_COVARIATE")),
        extractTraits.out.treatments.ifEmpty(file("NO_EXTRA_TREATMENT"))
    )

    // generate sieve estimates
    generateSieveEstimates(generateEstimates.out, generateIIDGenotypes.out)

    // generate Summaries
    generateSummaries(generateEstimates.out, generateSieveEstimates.out)
}