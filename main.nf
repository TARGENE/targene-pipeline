#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.DECRYPTED_DATASET = "NO_FILE"
params.CALL_THRESHOLD = 0.9
params.POSITIVITY_CONSTRAINT = 0.01
params.MAF_THRESHOLD = 0.01
params.VERBOSITY = 1
params.COHORT = "UKBB"
params.TRAITS_CONFIG = "NO_UKB_TRAIT_CONFIG"
params.WITHDRAWAL_LIST = 'NO_WITHDRAWAL_LIST'
params.OUTDIR = "${launchDir}/results"

// Confounding adjustment by PCA
params.NB_PCS = 6
params.QC_FILE = "NO_QC_FILE"
params.LD_BLOCKS = "NO_LD_BLOCKS"
params.FLASHPCA_EXCLUSION_REGIONS = "data/exclusion_regions_hg19.txt"

// Estimands Generation 
params.BATCH_SIZE = 400

// CUSTOM
params.STUDY_DESIGN = "CUSTOM"
params.ESTIMANDS_FILE = "NO_ESTIMANDS_FILE"

// FROM_ACTORS
params.EXTRA_CONFOUNDERS = 'NO_EXTRA_CONFOUNDER'
params.EXTRA_COVARIATES = 'NO_EXTRA_COVARIATE'
params.ENVIRONMENTALS = 'NO_EXTRA_TREATMENT'
params.ORDERS = "1,2"

// Sieve Variance Plateau
params.SVP = false
params.GRM_NSPLITS = 100
params.NB_SVP_ESTIMATORS = 100
params.MAX_SVP_THRESHOLD = 0.9
params.SVP_ESTIMATOR_KEY = "TMLE"

// TMLE
params.KEEP_IC = params.SVP == true ? true : false
params.PVAL_THRESHOLD = 0.05
params.TMLE_SAVE_EVERY = 100

// Outputs
params.ARROW_OUTPUT = "dataset.arrow"
params.JSON_OUTPUT = "NO_JSON_OUTPUT"
params.HDF5_OUTPUT = "results.hdf5"

// Permutation Tests Parameters
params.MAX_PERMUTATION_TESTS = null
params.PVAL_COL = "TMLE_PVALUE"
params.PERMUTATION_ORDERS = "1"
params.RNG = 123
params.MAF_MATCHING_RELTOL = 0.05
params.N_RANDOM_VARIANTS = 10


include { IIDGenotypes } from './modules/genotypes.nf'
include { FlashPCA; AdaptFlashPCA } from './modules/confounders.nf'
include { UKBFieldsList; UKBConv; TraitsFromUKB } from './modules/ukb_traits.nf'
include { TMLE; TMLEInputsFromParamFile; TMLEInputsFromActors } from './modules/tmle.nf'
include { GRMPart; AggregateGRM } from './modules/grm.nf'
include { SieveVarianceEstimation ; MergeOutputs } from './modules/sieve_variance.nf'
include { GeneratePermutationTestsData; GenerateRandomVariantsTestsData } from './modules/negative_control.nf'

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

    if (params.COHORT == "UKBB") {
        extracted_traits = TraitsFromUKB(decrypted_dataset, traits_config, withdrawal_list)
    } 
    else {
        extracted_traits = decrypted_dataset 
    }

    emit:
        extracted_traits
}

workflow generateIIDGenotypes {
    take:
        traits

    main:
        qc_file = Channel.value(file("$params.QC_FILE"))
        flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
        ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
        bed_files_ch = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

        IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file, traits)

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

workflow runPCA {
    // Extract traits
    extractTraits()

    // Generate IID Genotypes
    generateIIDGenotypes(extractTraits.out)

    // Genetic confounders up to NB_PCS
    geneticConfounders(generateIIDGenotypes.out) 

}

workflow generateTMLEEstimates {
    take:
        traits
        genetic_confounders

    main:
        estimator_file = Channel.value(file("$params.ESTIMATOR_FILE", checkIfExists: true))
        bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()

        if (params.STUDY_DESIGN == "FROM_ACTORS") {
            bqtls = Channel.value(file("$params.BQTLS"))
            trans_actors = Channel.fromPath("$params.TRANS_ACTORS", checkIfExists: true).collect()
            extra_confounders = Channel.value(file("$params.EXTRA_CONFOUNDERS"))
            extra_treatments = Channel.value(file("$params.ENVIRONMENTALS"))
            extra_covariates = Channel.value(file("$params.EXTRA_COVARIATES"))
            tmle_inputs = TMLEInputsFromActors(
                bgen_files,
                traits,
                genetic_confounders,
                extra_confounders,
                extra_treatments,
                extra_covariates,
                bqtls,
                trans_actors)
        }
        else if (params.STUDY_DESIGN == "CUSTOM"){
            estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                estimands_file)
        }
        else { 
            throw new Exception("This STUDY_DESIGN is not available.")
        }
        // compute TMLE estimates for continuous targets
        TMLE(
            tmle_inputs.traits,
            tmle_inputs.estimands.flatten(),
            estimator_file,
        )
        

    emit:
        TMLE.out
}


workflow generateSieveEstimates {
    take:
        tmle_files
        iid_genotypes
    
    main:
        grm_parts = Channel.from( 1..params.GRM_NSPLITS )
        GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)
        AggregateGRM(GRMPart.out.collect())
        // Sieve estimation
        SieveVarianceEstimation(tmle_files.collect(), AggregateGRM.out.grm_ids, AggregateGRM.out.grm_matrix)
    emit:
        SieveVarianceEstimation.out
}

workflow negativeControl {
    results_file = Channel.value(file("${params.HDF5_OUTPUT}"))

    // Permutation Tests
    dataset = Channel.value(file("${params.OUTDIR}/${params.ARROW_OUTPUT}"))
    estimator_file = Channel.value(file("${params.ESTIMATOR_FILE}"))
    sieve_csv = Channel.value(file("NO_SIEVE_FILE"))
    GeneratePermutationTestsData(dataset, results_file)
    TMLE(
        GeneratePermutationTestsData.output.dataset, 
        GeneratePermutationTestsData.output.parameters.flatten(), 
        estimator_file
    )
    MergeOutputs(TMLE.out.tmle_csv.collect(), sieve_csv, "permutation_summary.csv")
    
    // Random Variants parameter files generation
    if (params.STUDY_DESIGN == "FROM_ACTORS") {
        bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS", checkIfExists: true).collect()
        GenerateRandomVariantsTestsData(trans_actors, bgen_files, results_file)
    }
}

workflow {
    // Extract traits for UKBB
    extractTraits()
    
    // Generate IID Genotypes
    generateIIDGenotypes(extractTraits.out)

    // Genetic confounders
    geneticConfounders(generateIIDGenotypes.out)

    // generate estimates
    generateTMLEEstimates(
        extractTraits.out,
        geneticConfounders.out,
    )

    // Generate sieve estimates
    if (params.SVP == true){
        sieve_results = generateSieveEstimates(generateTMLEEstimates.out, generateIIDGenotypes.out)
    }

    MergeOutputs(generateTMLEEstimates.out.collect(), "summary.csv")
}
