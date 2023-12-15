#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.VERBOSITY = 0 
params.TRAITS_DATASET = "You need to provide a Traits dataset."
params.UKB_ENCODING_FILE = "NO_UKB_ENCODING_FILE"
params.CALL_THRESHOLD = 0.9
params.POSITIVITY_CONSTRAINT = 0.01
params.MAF_THRESHOLD = 0.01
params.COHORT = "UKBB"
params.UKB_CONFIG = "data/ukbconfig.yaml"
params.UKB_WITHDRAWAL_LIST = 'data/NO_WITHDRAWAL_LIST'
params.OUTDIR = "${launchDir}/results"

// Confounding adjustment by PCA
params.NB_PCS = 6
params.QC_FILE = "data/NO_QC_FILE"
params.LD_BLOCKS = "data/NO_LD_BLOCKS"
params.FLASHPCA_EXCLUSION_REGIONS = "data/exclusion_regions_hg19.txt"

// Estimands Generation 
params.BATCH_SIZE = 400

// CUSTOM
params.STUDY_DESIGN = "CUSTOM"
params.ESTIMANDS_FILE = "NO_ESTIMANDS_FILE"

// FROM_ACTORS
params.EXTRA_CONFOUNDERS = 'data/NO_EXTRA_CONFOUNDER'
params.EXTRA_COVARIATES = 'data/NO_EXTRA_COVARIATE'
params.ENVIRONMENTALS = 'data/NO_EXTRA_TREATMENT'
params.ORDERS = "1,2"

// SVP Parameters
params.SVP = false
params.GRM_NSPLITS = 100
params.NB_SVP_ESTIMATORS = 100
params.MAX_SVP_THRESHOLD = 0.9
params.SVP_ESTIMATOR_KEY = "TMLE"

// TMLE Parameters
params.KEEP_IC = params.SVP == true ? true : false
params.PVAL_THRESHOLD = 0.05
params.TMLE_SAVE_EVERY = 100
params.AGGREGATED_DATASET = "results/dataset.arrow"
params.ESTIMATOR_FILE = "glmnet"

// Outputs Parameters
params.ARROW_OUTPUT = "dataset.arrow"
params.JSON_OUTPUT = "NO_JSON_OUTPUT"
params.HDF5_OUTPUT = "results.hdf5"

include { EstimationInputs } from './modules/estimation_inputs.nf'
include { IIDGenotypes } from './modules/genotypes.nf'
include { GeneticConfounders } from './modules/confounders.nf'
include { ExtractTraits } from './modules/traits.nf'
include { EstimationWorkflow } from './modules/estimation.nf'
include { SVPWorkflow } from './modules/svp.nf'

log.info """\
         ${workflow.manifest.name} v${workflow.manifest.version}
         ==========================
         Cohort Type  : ${params.COHORT}
         Study Design : ${params.STUDY_DESIGN}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}
         """
         .stripIndent()


workflow {
    // Define Parameters
    verbosity = params.VERBOSITY

    study_design = params.STUDY_DESIGN
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    bqtls_file = Channel.value(file("$params.BQTLS"))
    transactors_files = Channel.fromPath("$params.TRANS_ACTORS").collect()
    extra_confounders = Channel.value(file("$params.EXTRA_CONFOUNDERS"))
    extra_treatments = Channel.value(file("$params.ENVIRONMENTALS"))
    extra_covariates = Channel.value(file("$params.EXTRA_COVARIATES"))
    batch_size = params.BATCH_SIZE
    call_threshold = params.CALL_THRESHOLD
    positivity_constraint = params.POSITIVITY_CONSTRAINT

    cohort = params.COHORT
    ukb_encoding_file = params.UKB_ENCODING_FILE
    ukb_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

    maf_threshold = params.MAF_THRESHOLD
    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    
    estimator_config = Channel.value(file("${params.ESTIMATOR_FILE}"))
    keep_ic = params.KEEP_IC
    pval_threshold = params.PVAL_THRESHOLD
    save_every = params.TMLE_SAVE_EVERY
    hdf5_output = "${params.HDF5_OUTPUT}"
    json_output = "${params.JSON_OUTPUT}"

    do_svp = params.SVP
    n_svp_estimators = params.NB_SVP_ESTIMATORS
    max_svp_threshold = params.MAX_SVP_THRESHOLD
    svp_estimator_key = params.SVP_ESTIMATOR_KEY
    grm_n_splits = params.GRM_NSPLITS

    // Extract Traits
    ExtractTraits(
        traits_dataset,
        cohort,
        ukb_config,
        ukb_withdrawal_list,
        ukb_encoding_file,
    )
    
    // IID Genotypes
    IIDGenotypes(
        flashpca_excl_reg,
        ld_blocks,
        bed_files,
        qc_file,
        ExtractTraits.out,
        maf_threshold
    )

    // Genetic confounders
    GeneticConfounders(IIDGenotypes.out)

    // generate main dataset and estimand configuration files
    EstimationInputs(
        study_design,
        bgen_files,
        ExtractTraits.out,
        GeneticConfounders.out,
        estimands_file,
        bqtls_file,
        transactors_files,
        extra_confounders,
        extra_treatments,
        extra_covariates,
        batch_size,
        call_threshold,
        positivity_constraint
    )

    // generate estimates
    EstimationWorkflow(
        EstimationInputs.out.aggregated_dataset,
        EstimationInputs.out.estimands.flatten(),
        estimator_config,
        keep_ic,
        do_svp,
        pval_threshold,
        save_every,
        hdf5_output,
        json_output
    )

    // Generate sieve estimates
    if (do_svp == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result, 
            IIDGenotypes.out,
            n_svp_estimators,
            max_svp_threshold,
            svp_estimator_key,
            grm_n_splits,
            verbosity
        )
    }
}
