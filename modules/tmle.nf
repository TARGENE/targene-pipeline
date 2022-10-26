def longest_prefix(files){
    // Only one file, strangely it is not passed as a list
    if (files instanceof Collection == false) {
        return files.getName()
    }
    // More than one file
    index = 0
    while(true){
        current_prefix = files[0].getName()[0..index]
        for (file in files){
            if(file.getName()[0..index] != current_prefix){
                return current_prefix[0..-2]
            }
        }
        index++
    }
}

process TMLE {
    container "olivierlabayle/targeted-estimation:new_interface"
    publishDir "$params.OUTDIR/summaries",  mode: 'symlink', pattern: "*.csv"
    publishDir "$params.OUTDIR/hdf5files",  mode: 'symlink', pattern: "*.hdf5"
    label "bigmem"
    label "multithreaded"

    input:
        path data
        path parameterfile
        path estimatorfile
    
    output:
        path "$outfilename"
    
    script:
        save_full = params.NB_VAR_ESTIMATORS !== 0 ? '--save-full' : ''
        outfilename = params.NB_VAR_ESTIMATORS !== 0 ? parameterfile.getName().replace("yaml", "hdf5") : parameterfile.getName().replace("yaml", "csv")
        """
        julia --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/scripts/tmle.jl \
        $data $parameterfile $estimatorfile $outfilename $save_full
        """
}

process TMLEInputsFromParamFiles {
    container "olivierlabayle/tl-core:new_strategies"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path parameters

    output:
        path "final.data.csv", emit: traits
        path "final.*.yaml", emit: parameters

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        params_prefix = longest_prefix(parameters)
        batch_size = params.PHENOTYPES_BATCH_SIZE == 0 ? "" :  "--phenotype-batch-size ${params.PHENOTYPES_BATCH_SIZE}"
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        $batch_size \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        from-param-files $params_prefix
        """
}

process TMLEInputsFromActors {
    container "olivierlabayle/tl-core:new_strategies"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path extra_confounders
        path extra_treatments
        path extra_covariates
        path bqtls
        path trans_actors

    output:
        path "final.data.csv", emit: traits
        path "final.*.yaml", emit: parameters

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        trans_actors_prefix = longest_prefix(trans_actors)
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        extra_covariates = extra_covariates.name != 'NO_EXTRA_COVARIATE' ? "--extra-covariates $covariates" : ''
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        --phenotype-batch-size ${params.PHENOTYPES_BATCH_SIZE} \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        from-actors $bqtls $trans_actors_prefix $extra_confounders $extra_treatments $extra_covariates --orders ${params.ORDERS}
        """
}
