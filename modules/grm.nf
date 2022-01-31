process GRMPart {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    label "bigmem"

    input:
        path bedfile
        val nparts
        val part_id

    output:
        path "UKKK.*.grm.*"
    
    script:
        "gcta64 --bfile $bedfile --make-grm-part $nparts $part_id --thread-num ${task.cpus} --out UKKK"
}

process AggregateGRMFiles {
    input:
        path ukbb_grm_files
        val extension

    output:
        path "UKBB$extension"

    script:
        "cat UKBB*$extension > UKBB$extension"
}