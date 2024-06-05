include { PCA } from './pca.nf'
include { EstimationInputs } from '../subworkflows/estimation_inputs.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'

workflow TARGENE {
    // Define Parameters
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    estimator_config = Channel.value(file("${params.ESTIMATOR_FILE}"))

    // DEPRECATED: LEGACY TO BE REMOVED
    bqtls_file = Channel.value(file("$params.BQTLS"))
    transactors_files = Channel.fromPath("$params.TRANS_ACTORS").collect()
    extra_confounders = Channel.value(file("$params.EXTRA_CONFOUNDERS"))
    extra_treatments = Channel.value(file("$params.ENVIRONMENTALS"))
    extra_covariates = Channel.value(file("$params.EXTRA_COVARIATES"))

    // Extract Traits
    PCA()

    // generate main dataset and estimand configuration files
    EstimationInputs(
        bgen_files,
        PCA.out.traits,
        PCA.out.pcs,
        estimands_file,
        bqtls_file,
        transactors_files,
        extra_confounders,
        extra_treatments,
        extra_covariates,
    )

    // generate estimates
    EstimationWorkflow(
        EstimationInputs.out.aggregated_dataset,
        EstimationInputs.out.estimands.flatten(),
        estimator_config,
    )

    // Generate sieve variance plateau estimates
    if (params.SVP == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result, 
            PCA.out.iid_genotypes,
        )
    }
}