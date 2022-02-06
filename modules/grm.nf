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

process AggregateGRMFiles {
    publishDir "$params.OUTDIR/GRM", mode: 'symlink'

    input:
        path grmfiles

    output:
        path "UKBB.Full*"

    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/aggregate_grm.jl UKBB UKBB_GRM"
}