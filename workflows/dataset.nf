include { MakeDataset } from '../modules/dataset.nf'
include { IIDGenotypes } from '../modules/genotypes.nf'
include { GeneticConfounders } from '../modules/confounders.nf'
include { ExtractTraits } from '../modules/traits.nf'

workflow MAKE_DATASET {
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
    bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()

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

    // Dataset
    MakeDataset(
        bgen_files,
        ExtractTraits.out,
        GeneticConfounders.out,
    )
}