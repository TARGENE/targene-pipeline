include { PCA } from './pca.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'

workflow TARGENE {
    // Define Parameters
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    estimator_config = Channel.value(file("${params.ESTIMATOR_FILE}"))

    // Extract Traits
    PCA()

    // generate main dataset and estimand configuration files
    pcs_and_genotypes = PCA.out.pcs.combine(bgen_files).map{it -> [it[0][0], it[0][1], it[1]]}
    EstimationInputs(
        pcs_and_genotypes,
        PCA.out.traits,
        estimands_file
    )

    // generate estimates
    EstimationWorkflow(
        EstimationInputs.out.dataset,
        EstimationInputs.out.estimands.flatten(),
        estimator_config,
    )

    // Generate sieve variance plateau estimates
    genotypes = PCA.out.iid_genotypes.map{it -> it[1]}
    if (params.SVP == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result, 
            genotypes,
        )
    }
}