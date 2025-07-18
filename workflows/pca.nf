include { IIDGenotypes; LOCOGenotypes; FlashPCA } from '../subworkflows/confounders.nf'
include { ExtractTraits } from '../subworkflows/extract_traits.nf'

workflow PCA {
    // Define Parameters
    ukb_encoding_file = params.UKB_ENCODING_FILE
    ukb_config = Channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
    ukb_withdrawal_list = Channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
    traits_dataset = Channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

    qc_file = Channel.value(file("$params.QC_FILE", checkIfExists: true))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
    bed_files = Channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    
    // Extract Traits
    ExtractTraits(
        traits_dataset,
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
    )

    // PCA
    FlashPCA(IIDGenotypes.out)
    
    emit:
        traits = ExtractTraits.out
        iid_genotypes = IIDGenotypes.out
        pcs = FlashPCA.out.pcs

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
    
    // Extract Traits
    ExtractTraits(
        traits_dataset,
        ukb_config,
        ukb_withdrawal_list,
        ukb_encoding_file,
    )

    // IID Genotypes
    loco_genotypes = LOCOGenotypes(
        flashpca_excl_reg,
        ld_blocks,
        bed_files, 
        qc_file,
        ExtractTraits.out
    )

    // PCA
    FlashPCA(loco_genotypes)

    emit:
        traits = ExtractTraits.out
        iid_genotypes = loco_genotypes
        confounders = FlashPCA.out.pcs
}
