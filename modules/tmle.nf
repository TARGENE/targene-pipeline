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
    publishDir "$params.OUTDIR/hdf5files", saveAs: { filename -> filename.split("_batch")[0] + "/$filename" }, mode: 'symlink'
    label "bigmem"
    label "multithreaded"

    input:
        path treatments
        path targets 
        path confounders
        path parameterfile
        path estimator
        path covariatesfile
        val target_type
    
    output:
        path "*.hdf5"
    
    script:
        covariates = covariatesfile.name != 'NO_COVARIATE' ? "--covariates $covariatesfile" : ''
        save_models = params.SAVE_MODELS == true ? '--save-models' : ''
        save_ic = params.SAVE_IC == true ? '' : '--no-ic'
        outfilename = parameterfile.getName().replace("yaml", "hdf5")
        """
        julia --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/scripts/tmle.jl \
        $treatments $targets $confounders $parameterfile $estimator $outfilename \
        $covariates --target-type $target_type $save_models $save_ic
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
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        --phenotype-batch-size ${params.PHENOTYPES_BATCH_SIZE} \
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
