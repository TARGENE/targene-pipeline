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

process TMLEInputsFromGivenParams {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"

    input:
        path bgenfiles
        path binary_phenotypes
        path continuous_phenotypes
        path genetic_confounders
        path extra_confounders
        path extra_treatments
        path covariates
        path parameters

    output:
        path "final.binary-phenotypes.csv", emit: binary_phenotypes
        path "final.continuous-phenotypes.csv", emit: continuous_phenotypes
        path "final.confounders.csv", emit: confounders
        path "final.covariates.csv", emit: covariates, optional: true
        path "final.treatments.csv", emit: treatments
        path "final.binary.*.yaml", emit: binary_parameters
        path "final.continuous.*.yaml", emit: continuous_parameters

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        params_prefix = longest_prefix(parameters)
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        covariates = covariates.name != 'NO_COVARIATE' ? "--covariates $covariates" : ''
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --binary-phenotypes $binary_phenotypes --continuous-phenotypes $continuous_phenotypes \
        --bgen-prefix $bgen_prefix --call-threshold ${params.CALL_THRESHOLD} \
        --genetic-confounders $genetic_confounders $extra_confounders \
        $extra_treatments $covariates \
        --phenotype-batch-size ${params.PHENOTYPES_BATCH_SIZE} \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        with-param-files $params_prefix
        """
}

process TMLEInputsFromASBTrans {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"

    input:
        path bgenfiles
        path binary_phenotypes
        path continuous_phenotypes
        path genetic_confounders
        path extra_confounders
        path extra_treatments
        path covariates
        path asbs
        path trans_actors
        path parameters

    output:
        path "final.binary-phenotypes.csv", emit: binary_phenotypes
        path "final.continuous-phenotypes.csv", emit: continuous_phenotypes
        path "final.confounders.csv", emit: confounders
        path "final.covariates.csv", emit: covariates, optional: true
        path "final.treatments.csv", emit: treatments
        path "final.binary.*.yaml", emit: binary_parameters
        path "final.continuous.*.yaml", emit: continuous_parameters

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        asb_prefix = longest_prefix(asbs)
        param_prefix = parameters[0].name ? "NO_PARAMETER_FILE" : "--param-prefix ${longest_prefix(parameters)}"
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        covariates = covariates.name != 'NO_COVARIATE' ? "--covariates $covariates" : ''
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --binary-phenotypes $binary_phenotypes --continuous-phenotypes $continuous_phenotypes \
        --bgen-prefix $bgen_prefix --call-threshold ${params.CALL_THRESHOLD} \
        --genetic-confounders $genetic_confounders $extra_confounders \
        $extra_treatments $covariates \
        --phenotype-batch-size ${params.PHENOTYPES_BATCH_SIZE} \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        with-asb-trans $asb_prefix $trans_actors $param_prefix
        """
}
