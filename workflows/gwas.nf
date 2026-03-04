include { CreateEstimatorsConfigChannel } from '../modules/utils.nf'
include { LocoPCA } from './pca.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'
include { IIDGenotypes } from '../subworkflows/confounders.nf'

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
        EstimationInputs.out.transpose(), estimator_config
    )

    // Generate sieve variance plateau estimates
    if (params.SVP == true) {
        if (params.PREVALENCE != "NO_SET_PREVALENCE") {
            error "SVP is not compatible with a set PREVALENCE parameter."
        }
        // IID Genotypes for SVP (needs all chromosomes merged)
        qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
        flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
        ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
        IIDGenotypes(
            flashpca_excl_reg,
            ld_blocks,
            bed_files,
            qc_file,
            LocoPCA.out.traits,
        )
        genotypes = IIDGenotypes.out.map{genotypes_id, genotypes -> genotypes}.collect()
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result.collect(), 
            genotypes,
        )
    }
}