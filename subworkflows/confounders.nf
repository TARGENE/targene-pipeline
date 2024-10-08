include { filterBED; thinByLD; mergeBEDS; SampleQCFilter; FlashPCA } from '../modules/confounders.nf'
include { leave_chr_out } from '../modules/utils.nf'

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
            .map{it -> ["all_genotypes", it]}
        mergeBEDS(bedfiles_to_be_merged)
        SampleQCFilter(mergeBEDS.out.collect())

    emit:
        SampleQCFilter.out
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