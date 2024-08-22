include { filterBED; thinByLD; mergeBEDS; SampleQCFilter; FlashPCA; AdaptFlashPCA } from '../modules/confounders.nf'

workflow IIDGenotypes{
    take:
        flashpca_excl_reg
        ld_blocks
        bed_files
        qc_file
        traits

    main:
        filtered_bedfiles = filterBED(bed_files, qc_file, ld_blocks, traits)
        ld_pruned = thinByLD(flashpca_excl_reg, filtered_bedfiles)
        bedfiles_to_be_merged = ld_pruned.collect()
            .map{it -> ["ukbb_merged", it]}
        mergeBEDS(bedfiles_to_be_merged)
        SampleQCFilter(mergeBEDS.out.collect())

    emit:
        SampleQCFilter.out
}

workflow GeneticConfounders {
    take:
        iid_genotypes

    main:
        FlashPCA(iid_genotypes)
        AdaptFlashPCA(FlashPCA.out)

    emit:
        AdaptFlashPCA.out
}