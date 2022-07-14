// Would be good if we could seprate a dataframe of eQTLs/bQTLs with rsid, chr and position into seprate rows that maintain their dataframe format.
process  gcta_ld{
    container "olivierlabayle/ukbmain:v0.1.0"

    input:
        path snps_list
        path bed_files_ch

    output:
        path snps_in_ld

    script:
        """
        gcta64 --bfile data/ukb_53116_whole_genome --ld snp_list.txt --ld-wind 10000 --ld-sig 0.05 --out snps_in_ld -output
        """
}

process  position_gcta_snps{
    container "olivierlabayle/ukbmain:v0.1.0"

    input:
        path ld_list_of_snps.snp.ld
        path bed_files_ch

    output:
        path GCTA_position

    script:
        """
        plink --bfile ukb_53116_chr --extract snps_in_ld.snp.ld  --make-bed  --out GCTA_position

        """
}

process create_blocks{
    container "olivierlabayle/ukbmain:v0.1.0"
// I need to read in the file names as strings in the python script. Do you know how to do this with nextflow?
    input:
        path snps_in_ld
        path bed_files_ch
        path GCTA_position

    output:
        path "ld_position_of_snp.csv"

    script:
        """
        python Creating_LD_block.py

        """
}


workflow LD_blocks{
    list_snps = Channel.value(file("$params.snps_list"))
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    LD_blocks(list_snps, bed_files_ch)

    emit:
        ld_blocks.out
}