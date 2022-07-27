process TMLE {
    container "olivierlabayle/targeted-estimation:0.1.0"
    publishDir "$params.OUTDIR/hdf5files", saveAs: { filename -> filename.split("_batch")[0] + "/$filename" }, mode: 'symlink'
    label "bigmem"
    label "multithreaded"

    input:
        path genotypefile
        path phenotypefile
        path covariatesfile
        path estimatorfile
        tuple file(queryfile), file(phenotypes_batch)
        val target_type
    
    output:
        path "*.hdf5"
    
    script:
        save_full = params.SAVE_FULL == true ? '--save-full' : ''
        queryfilename = queryfile.getName()
        phen_batch = phenotypes_batch.getName()
        batch_id = phen_batch[17..-5]
        """
        outfilename=\$(julia --project --startup-file=no -e 'using TOML; ks=join(sort(collect(keys(TOML.parse(open("${queryfilename}"))["SNPS"]))), "_");println(ks)')
        outfilename="\${outfilename}_batch_${batch_id}_${target_type}.hdf5"
        julia --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/scripts/tmle.jl \
        $genotypefile $phenotypefile $covariatesfile $queryfile $estimatorfile \$outfilename \
        --phenotypes-list=$phen_batch --target-type=$target_type $save_full
        """
}


process PhenotypesBatches {
    container "olivierlabayle/tl-core:sample_filtering"

    input:
        path phenotypes_file
    
    output:
        path "phenotypes_batch_*"
    
    script:
        "julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/make_phenotypes_batches.jl ${phenotypes_file.getName()} --batch-size $params.PHENOTYPES_BATCH_SIZE"
}

process FinalizeTMLEInputs {
    input:
        path continuous_phenotypes
        path binary_phenotypes
        path genotypes
        path genetic_confounders
        path extra_confounders
        path covariates
        path extra_treatments
    
    output:
        path "final.binary-phenotypes.csv", emit: binary_phenotypes
        path "final.continuous-phenotypes.csv", emit: continuous_phenotypes
        path "final.treatments.csv", emit: treatments
        path "final.confounders.csv", emit: confounders
        path "final.covariates.csv", emit: covariates, optional: true
    
    script:
        covariates_option = covariates == "NO_FILE" ? "--covariates $covariates" : ""
        treatments_option = extra_treatments == "NO_FILE" ? "--extra_treatments $extra_treatments" : ""
        extra_confounders_option = extra_confounders == "NO_FILE" ? "--extra-confounders $extra_confounders" : ""
        """
        julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/finalize_tmle_inputs.jl \
        --out-prefix final \
        --binary-phenotypes $binary_phenotypes \
        --continuous-phenotypes $continuous_phenotypes \
        --genotypes $genotypes \
        --genetic-confounders $genetic_confounders \
        $extra_confounders_option $covariates_option $treatments_option
        """
}