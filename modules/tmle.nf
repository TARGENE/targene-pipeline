process TMLE {
    container "olivierlabayle/tmle-epistasis:0.3.1"
    publishDir "$params.OUTDIR/hdf5files", saveAs: { filename -> filename.split("_batch")[0] + "/$filename" }, mode: 'symlink'
    label "bigmem"
    label "multithreaded"

    input:
        path genotypefile
        path phenotypefile
        path confoundersfile
        path estimatorfile
        tuple file(queryfile), file(phenotypes_batch)
        val target_type
    
    output:
        path "*.hdf5"
    
    script:
        adaptive_cv = params.ADAPTIVE_CV == true ? '--adaptive-cv' : ''
        save_full = params.SAVE_FULL == true ? '--save-full' : ''
        queryfilename = queryfile.getName()
        phen_batch = phenotypes_batch.getName()
        batch_id = phen_batch[17..-5]
        """
        outfilename=\$(julia --project --startup-file=no -e 'using TOML; ks=join(sort(collect(keys(TOML.parse(open("${queryfilename}"))["SNPS"]))), "_");println(ks)')
        outfilename="\${outfilename}_batch_${batch_id}_${target_type}.hdf5"
        julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/ukbb.jl \
        $genotypefile $phenotypefile $confoundersfile $queryfile $estimatorfile \$outfilename \
        --phenotypes-list=$phen_batch --min-freq=${params.MINOR_CAT_FREQUENCY} --target-type=$target_type $adaptive_cv $save_full
        """
}


process PhenotypesBatches {
    container "olivierlabayle/tmle-epistasis:0.3.1"

    input:
        path phenotypes_file
    
    output:
        path "phenotypes_batch_*"
    
    script:
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/make_phenotypes_batches.jl ${phenotypes_file.getName()} --batch-size $params.PHENOTYPES_BATCH_SIZE"
}