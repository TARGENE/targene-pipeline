process GRM {
    container "olivierlabayle/ukbb-estimation-pipeline:0.2.1"
    label "verybigmem"

    input:
        path bedfiles
    
    output:
        path "grmmatrix"
    
    script:
        only_bedfiles = bedfiles.findAll { it.getName().endsWith(".bed") }
        "julia --threads=${task.cpus} --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/compute_grm.jl grmmatrix ${only_bedfiles.join(" ")}"
}