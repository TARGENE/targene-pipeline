process GRMPart {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    label "bigmem"
    label "multithreaded"

    input:
        path bedfiles
        val nparts
        val part_id

    output:
        path "UKBB*.grm.*"
    
    script:
        base = bedfiles.first().getName().split("\\.")[0]
        "gcta64 --bfile $base --make-grm-part $nparts $part_id --thread-num ${task.cpus} --out UKBB"
}

process AggregateIDFiles {
    publishDir "$params.OUTDIR/GRM", mode: 'symlink'

    input:
        path grm_id_files

    output:
        path "UKBB_GRM.ids"

    script:
        filenames = grm_id_files.collect { it.getName() }.join(" ")
        "cat $filenames > UKBB_GRM.ids"
}


process PrependSize {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"
    publishDir "$params.OUTDIR/GRM", mode: 'symlink'
    label "bigmem"

    input:
        path grmpart_file

    output:
        path "GRM_*"

    script:
        filename = grmpart_file.getName()
        file_number = filename.split("_")[-1][0..-9]
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/prepend_grmpart_size.jl $filename GRM_${file_number}.bin"
}