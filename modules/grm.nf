process GRM {
    container "olivierlabayle/ukbb-estimation-pipeline:0.2.1"
    label "verybigmem"

    input:
        path bedfiles
    
    output:
        path "grm.hdf5"
    
    script:
        only_bedfiles = bedfiles.findAll { it.getName().endsWith(".bed") }
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/compute_grm.jl grm.hdf5 ${only_bedfiles.join(" ")}"
}