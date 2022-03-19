def String build_outfilename(queryfile, phenotypes_batch, type) {
    allLines = queryfile.readLines()
    SNPSection = false
    outfilename = ""
    for (line : allLines) {
        if (line == "[SNPS]") {
            SNPSection = true
        }
        else if (SNPSection) {
            if (line.contains("=")) {
                newRS = line.split("=")[0].strip()
                outfilename += newRS + "_"
            }
        else if (line.startsWith("["))
            SNPSection = false
        }
    }

    outfilename += "batch_" +phenotypes_batch.getName()[17..-5] + "_" + type + ".hdf5"

    return outfilename
}

process TMLE {
    container "olivierlabayle/tmle-epistasis:0.3.0"
    label "bigmem"
    label "multithreaded"

    input:
        path bgenfiles
        path phenotypefile
        path confoundersfile
        path estimatorfile
        tuple file(queryfile), file(phenotypes_batch)
        val target_type
    
    output:
        path "*.hdf5"
    
    script:
        def adaptive_cv = params.ADAPTIVE_CV == true ? '--adaptive-cv' : ''
        def save_full = params.SAVE_FULL == true ? '--save-full' : ''
        outfilename = build_outfilename(file(queryfile), phenotype_batch, type)
        """
        julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/ukbb.jl $phenotypefile $confoundersfile $queryfile $estimatorfile $outfilename --phenotypes-list ${phenotypes_batch.getName()} --target-type $target_type $adaptive_cv $save_full
        """
}


process PhenotypesBatches {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"

    input:
        path phenotypes_file
    
    output:
        "phenotypes_batch_*"
    
    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/make_phenotypes_batches.jl ${phenotypes_file.getName()} --batch-size $params.PHENOTYPES_BATCH_SIZE"
}

process MergeHDF5 {
    container "olivierlabayle/ukbb-estimation-pipeline:0.3.0"

    input:
        tuple val rs_ids , path results_files

    output:
        "${rs_ids}.hdf5"

    script:
        "julia --project=/EstimationPipeline.jl --startup-file=no /EstimationPipeline.jl/bin/merge_jld2.jl ${rs_ids}"
}