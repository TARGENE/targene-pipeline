process GRMPart {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    label "bigmem"
    label "multithreaded"

    input:
        path bedfiles
        val nparts
        val part_id

    output:
        path "GRM*.grm.*"
    
    script:
        base = bedfiles.first().getName().split("\\.")[0]
        "gcta64 --bfile $base --make-grm-part $nparts $part_id --thread-num ${task.cpus} --out GRM"
}

process AggregateGRM {
    publishDir "$params.OUTDIR/GRM", mode: 'symlink'

    input:
        path grm_files

    output:
        path "GRM.id", emit: grm_ids
        path "GRM.bin", emit: grm_matrix

    script:
        """
        cat *.grm.id > GRM.id
        cat *.grm.bin > GRM.bin
        """
}
