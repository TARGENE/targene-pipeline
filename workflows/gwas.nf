include { CreateEstimatorsConfigChannel } from '../modules/utils.nf'
include { LocoPCA } from './pca.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'

workflow GWAS {
    // Define Parameters
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    estimands_file = Channel.value(file("$params.ESTIMANDS_CONFIG"))
    estimator_config = CreateEstimatorsConfigChannel(params.ESTIMATORS_CONFIG)

    // Loco PCA
    LocoPCA()
    
    // Estimation Inputs
    pcs_and_genotypes = LocoPCA.out.confounders.join(bed_files, failOnDuplicate: true)
    EstimationInputs(
        pcs_and_genotypes,
        LocoPCA.out.traits,
        estimands_file
    )

    // Estimation
    EstimationWorkflow(
        EstimationInputs.out.transpose(),
        estimator_config,
    )
}