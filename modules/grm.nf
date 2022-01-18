process GRM {
    container "olivierlabayle/ukbb-estimation-pipeline:0.2.0"
    label "verybigmem"
    cpus 10

    input:
        path bedfiles
    
    output:
        path "grmmatrix"
    
    script:
        only_bedfiles = bedfiles.flatten().findAll { it.endsWith(".bed") }
        "julia --threads=${task.cpus} --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/compute_grm.jl grmmatrix ${only_bedfiles.join(" ")}"
}