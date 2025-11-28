process subsetBED{
    label 'bigmem'
    label 'plink_image'

    input:
        tuple val(chr_id), file(bedfiles)
        path subset_ids

    output:
        tuple val(chr_id), path("subset.${chr_id}.*"), emit: subset_bed_files

    script:
        input_prefix = bedfiles[0].toString().minus('.bed')
        """
        plink2 --bfile ${input_prefix} --extract ${subset_ids} --make-bed --out subset.${chr_id}
        """
}

process denseBED{
    label 'bigmem'
    label 'plink_image'

    input:
        tuple val(chr_id), file(bgenfiles)
        path prioritized_variants

    output:
        tuple val(chr_id), path("dense.${chr_id}.*"), emit: dense_bed_files

    script:
        input_prefix = bgenfiles[0].toString().minus('.bgen')
        """
        awk -v W=${params.WINDOW_SIZE} 'BEGIN { OFS = "\\t" }
            NR > 1 {
                start = \$2 - W
                if (start < 1) start = 1
                end = \$2 + W
                print \$1, start, end, \$3
            }' ${prioritized_variants} > ranges.${chr_id}.txt

        plink2 --bgen ${input_prefix}.bgen ref-first --sample ${input_prefix}.sample --make-bed --extract range ranges.${chr_id}.txt --out dense.${chr_id}
        """
}