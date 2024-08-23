
include { PCA } from './pca.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'
include { ExtractTraits } from '../subworkflows/extract_traits.nf'
include { filterBED; thinByLD; mergeBEDS; SampleQCFilter; FlashPCA; AdaptFlashPCA } from '../modules/confounders.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'

def filepath_matches_chr_prefix(fp, chr_prefix){
    def fp_string = fp.normalize().toString()
    return fp_string.contains(chr_prefix.normalize().toString() + ".")
}

def leave_chr_out(chr_prefix, bed_files){
    def bed_files_not_matching_chr_prefix = bed_files.findAll{ fp -> !filepath_matches_chr_prefix(fp, chr_prefix) }
    return [chr_prefix, bed_files_not_matching_chr_prefix]
}

workflow LOCOGenotypes {
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
        traits
        
    main:
        filtered_bed_files = filterBED(bed_files, qc_file, ld_blocks, traits)
        ld_pruned = thinByLD(flashpca_excl_reg, filtered_bed_files)
        bed_files_sets_to_be_merged = bed_files
            .combine(ld_pruned.collect().toList())
            .map {it -> leave_chr_out(it[0], it[2])}
        mergeBEDS(bed_files_sets_to_be_merged)
        SampleQCFilter(mergeBEDS.out)
    
    emit:
        SampleQCFilter.out
}

workflow LocoPCA {
    // Define Parameters
    ukb_encoding_file = params.UKB_ENCODING_FILE
    ukb_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    ExtractTraits(
        traits_dataset,
        ukb_config,
        ukb_withdrawal_list,
        ukb_encoding_file,
    )

    loco_genotypes = LOCOGenotypes(
        flashpca_excl_reg,
        ld_blocks,
        bed_files, 
        qc_file,
        ExtractTraits.out
    )
    // Genetic confounders
    FlashPCA(loco_genotypes)

    emit:
        traits = ExtractTraits.out
        iid_genotypes = loco_genotypes
        confounders = FlashPCA.out
}

workflow GWAS {
    // Define Parameters
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
    estimator_config = Channel.value(file("$params.ESTIMATOR_FILE"))

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