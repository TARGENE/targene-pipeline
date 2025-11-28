include { CreateEstimatorsConfigChannel } from '../modules/utils.nf'
include { LocoPCA } from './pca.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'
include { subsetBED; denseBED } from '../modules/extract_variants.nf'

workflow GWAS {
    // Define Parameters
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    estimands_file = Channel.value(file("$params.ESTIMANDS_CONFIG"))
    estimator_config = Channel.value(file("$params.ESTIMATORS_CONFIG"))

    // Loco PCA
    LocoPCA()
    
    // Estimation Inputs
    if (params.SUBSET_FILE != "NO_SUBSET_FILE") {
        subset_ids_file = Channel.value(file("$params.SUBSET_FILE", checkIfExists: true))
        target_bed_files = Channel.fromFilePairs("$params.TARGET_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
        subset_bed_files = subsetBED(target_bed_files, subset_ids_file).subset_bed_files
        
        pcs_and_genotypes = LocoPCA.out.confounders.join(subset_bed_files, failOnDuplicate: true)
    } else if (params.DENSE_MAPPING_FILE != "NO_DENSE_MAPPING_FILE") {
        prioritized_variants = Channel.value(file("$params.DENSE_MAPPING_FILE", checkIfExists: true))
        imputed_bgen_files = Channel.fromFilePairs("$params.BGEN_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
        dense_bed_files = denseBED(imputed_bgen_files, prioritized_variants).dense_bed_files

        pcs_and_genotypes = LocoPCA.out.confounders.join(dense_bed_files, failOnDuplicate: true)
    } else {
        pcs_and_genotypes = LocoPCA.out.confounders.join(bed_files, failOnDuplicate: true)
    }

    EstimationInputs(
        pcs_and_genotypes,
        LocoPCA.out.traits,
        estimands_file
    )

    // Estimation
    EstimationWorkflow(
        EstimationInputs.out.inputs.transpose(), estimator_config
    )
}