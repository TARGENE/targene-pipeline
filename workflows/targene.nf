include { PCA } from './pca.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'

workflow TARGENE {
    // Define Parameters
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect().toList()
    estimands_file = Channel.value(file("$params.ESTIMANDS_CONFIG"))
    estimator_config = Channel.value(file("${params.ESTIMATOR_CONFIG}"))

    // PCA
    PCA()

    // Estimation Inputs
    pcs_and_genotypes = PCA.out.pcs.combine(bgen_files)
    EstimationInputs(
        pcs_and_genotypes,
        PCA.out.traits,
        estimands_file
    )

    // generate estimates
    EstimationWorkflow(
        EstimationInputs.out.transpose(),
        estimator_config,
    )

    // Generate sieve variance plateau estimates
    genotypes = PCA.out.iid_genotypes.map{genotypes_id, genotypes -> genotypes}.collect()
    if (params.SVP == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result.collect(), 
            genotypes,
        )
    }
}