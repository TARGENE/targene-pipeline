#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.DECRYPTED_DATASET = "NO_FILE"
params.PARAMETER_PLAN = "FROM_PARAM_FILES"
params.PARAMETER_FILE = "NO_PARAMETER_FILE"
params.CALL_THRESHOLD = 0.9
params.POSITIVITY_CONSTRAINT = 0.01
params.SAVE_IC = true

params.MAF_THRESHOLD = 0.01
params.NB_PCS = 6

params.GRM_NSPLITS = 100
params.NB_VAR_ESTIMATORS = 0
params.MAX_TAU = 0.9
params.PVAL_THRESHOLD = 0.05

params.VERBOSITY = 1
params.OUTDIR = "${launchDir}/results"
params.RESULTS_FILE = "${params.OUTDIR}/summary.csv"
params.TMLE_INPUT_DATASET = "${params.OUTDIR}/tmle_inputs/final.data.csv"

params.TRAITS_CONFIG = "NO_UKB_TRAIT_CONFIG"
params.WITHDRAWAL_LIST = 'NO_WITHDRAWAL_LIST'

params.BATCH_SIZE = 400
params.EXTRA_CONFOUNDERS = 'NO_EXTRA_CONFOUNDER'
params.EXTRA_COVARIATES = 'NO_EXTRA_COVARIATE'
params.ENVIRONMENTALS = 'NO_EXTRA_TREATMENT'
params.ORDERS = "1,2"
params.GENOTYPES_AS_INT = false

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

    TraitsFromUKB(decrypted_dataset, traits_config, withdrawal_list)

    emit:
        TraitsFromUKB.out
}

workflow generateIIDGenotypes {
    take:
        traits

    main:
        qc_file = Channel.value(file("$params.QC_FILE"))
        flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
        ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
        bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

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

workflow generateTMLEEstimates {
    take:
        traits
        genetic_confounders

    main:
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))
        bgen_files = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true).collect()

        if (params.PARAMETER_PLAN == "FROM_ACTORS") {
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
        else if (params.PARAMETER_PLAN == "FROM_PARAM_FILE"){
            parameter_file = Channel.value(file("$params.PARAMETER_FILE"))
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                parameter_file)
        }
        else { 
            throw new Exception("This PARAMETER_PLAN is not available.")
        }
        // compute TMLE estimates for continuous targets
        TMLE(
            tmle_inputs.traits,
            tmle_inputs.parameters.flatten(),
            estimator_file,
        )
        

    emit:
        tmle_csvs = TMLE.out.tmle_csv
        inf_curves = TMLE.out.inf_curve
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
        csv_file = SieveVarianceEstimation.out.csv_file
        hdf5_file = SieveVarianceEstimation.out.hdf5_file
}

workflow negativeControl {
    results_file = Channel.value(file("${params.RESULTS_FILE}"))

    // Permutation Tests
    dataset = Channel.value(file("${params.TMLE_INPUT_DATASET}"))
    estimatorfile = Channel.value(file("${params.ESTIMATORFILE}"))
    sieve_csv = Channel.value(file("NO_SIEVE_FILE"))
    GeneratePermutationTestsData(dataset, results_file)
    TMLE(
        GeneratePermutationTestsData.output.dataset, 
        GeneratePermutationTestsData.output.parameters, 
        estimatorfile
    )
    MergeOutputs(TMLE.out.tmle_csv.collect(), sieve_csv, "permutation_summary.csv")
    
    // Random Variants parameter files generation
    if (params.PARAMETER_PLAN == "FROM_ACTORS") {
        bgen_files = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true).collect()
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS", checkIfExists: true).collect()
        GenerateRandomVariantsTestsData(trans_actors, bgen_files, results_file)
    }
}

workflow {
    // Extract traits
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

    // generate sieve estimates
    if (params.NB_VAR_ESTIMATORS != 0){
        sieve_results = generateSieveEstimates(generateTMLEEstimates.out.inf_curves, generateIIDGenotypes.out)
        sieve_csv = sieve_results.csv_file
    }
    else {
        sieve_csv = Channel.value(file("NO_SIEVE_FILE"))
    }

    MergeOutputs(generateTMLEEstimates.out.tmle_csvs.collect(), sieve_csv, "summary.csv")
}