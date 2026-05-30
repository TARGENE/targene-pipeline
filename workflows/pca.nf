include { IIDGenotypes; LOCOGenotypes } from '../subworkflows/confounders.nf'
include { FlashPCA; ProjectPCA } from '../modules/confounders.nf'
include { ExtractTraits } from '../subworkflows/extract_traits.nf'

// Add aliases to call the same process within the workflow
include { ExtractTraits as ExtractTraitsProjection } from '../subworkflows/extract_traits.nf'
include { IIDGenotypes as IIDGenotypesProjection } from '../subworkflows/confounders.nf'

workflow PCA {
    main:
        // Define Parameters
        ukb_encoding_file = params.UKB_ENCODING_FILE
        ukb_config = channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
        ukb_withdrawal_list = channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
        traits_dataset = channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))
        projection_dataset = channel.value(file("$params.PROJECTION_DATASET"))

        qc_file = channel.value(file("$params.QC_FILE", checkIfExists: true))
        flashpca_excl_reg = channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
        ld_blocks = channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
        bed_files = channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
        
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
            "all_genotypes",
        )

        // PCA
        FlashPCA(IIDGenotypes.out)

        // Optional Projection
        if (params.PROJECTION_DATASET != "NO_PROJECTION_DATASET") {
            // Extract Traits
            ExtractTraitsProjection(
                projection_dataset,
                ukb_config,
                ukb_withdrawal_list,
                ukb_encoding_file,
            )
            
            // IID Genotypes
            IIDGenotypesProjection(
                flashpca_excl_reg,
                ld_blocks,
                bed_files,
                qc_file,
                ExtractTraitsProjection.out,
                "projection_genotypes",
            )

            // PCA Projection
            projection_input = IIDGenotypesProjection.out
                .join(FlashPCA.out.loadings)
                .join(FlashPCA.out.meansd)
            ProjectPCA(projection_input)
            projected_pcs = ProjectPCA.out
        } else {
            projected_pcs = Channel.empty()
        }
    
    emit:
        traits = ExtractTraits.out
        iid_genotypes = IIDGenotypes.out
        pcs = FlashPCA.out.pcs
        projected_pcs = projected_pcs

}

workflow LocoPCA {
    main:
        // Define Parameters
        ukb_encoding_file = params.UKB_ENCODING_FILE
        ukb_config = channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
        ukb_withdrawal_list = channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
        traits_dataset = channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

        qc_file = channel.value(file("$params.QC_FILE", checkIfExists: true))
        flashpca_excl_reg = channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
        ld_blocks = channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
        bed_files = channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
        
        // Extract Traits
        ExtractTraits(
            traits_dataset,
            ukb_config,
            ukb_withdrawal_list,
            ukb_encoding_file,
        )

        // LOCO Genotypes for PCA (to avoid proximal contamination)
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
        confounders = FlashPCA.out.pcs
}
