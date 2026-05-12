include { LocoPCA } from './pca.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'
include { IIDGenotypes } from '../subworkflows/confounders.nf'
include { subsetBED; denseBED } from '../modules/extract_variants.nf'

workflow GWAS {
    // Define Parameters
    bed_files = channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    estimands_file = channel.value(file("$params.ESTIMANDS_CONFIG"))
    estimator_config = channel.fromPath(EstimatorsConfig.create(params.ESTIMATORS_CONFIG, params.OUTDIR))

    // Loco PCA
    LocoPCA()
    
    // Estimation Inputs
    if (params.SUBSET_FILE != "NO_SUBSET_FILE") {
        subset_ids_file = channel.value(file("$params.SUBSET_FILE", checkIfExists: true))
        target_bed_files = channel.fromFilePairs("$params.TARGET_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
        subset_bed_files = subsetBED(target_bed_files, subset_ids_file).subset_bed_files
        
        pcs_and_genotypes = LocoPCA.out.confounders.join(subset_bed_files, failOnDuplicate: true)
    } else if (params.DENSE_MAPPING_FILE != "NO_DENSE_MAPPING_FILE") {
        prioritized_variants = channel.value(file("$params.DENSE_MAPPING_FILE", checkIfExists: true))
        imputed_bgen_files = channel.fromFilePairs("$params.BGEN_FILES", size: 3, checkIfExists: true){ f -> f.name.replaceAll(/\.bgen(\.bgi)?$|\.sample$/,'') }
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
        EstimationInputs.out.transpose(), estimator_config
    )

    // Generate sieve variance plateau estimates
    if (params.SVP == true) {
        if (params.PREVALENCE != "NO_SET_PREVALENCE") {
            error "SVP is not compatible with a set PREVALENCE parameter."
        }
        // IID Genotypes for SVP (needs all chromosomes merged)
        qc_file = channel.value(file("$params.QC_FILE", checkIfExists: true))
        flashpca_excl_reg = channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
        ld_blocks = channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
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