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
    container "olivierlabayle/targeted-estimation:v0.3.0"
    publishDir "$params.OUTDIR/csvs",  mode: 'symlink', pattern: "*.csv"
    publishDir "$params.OUTDIR/hdf5files/inf_curves",  mode: 'symlink', pattern: "*.hdf5"
    label "bigmem"
    label "multithreaded"

    input:
        path data
        path parameterfile
        path estimatorfile
    
    output:
        path "${outprefix}.csv", emit: tmle_csv
        path "${outprefix}.hdf5", optional: true, emit: inf_curve
    
    script:
        save_ic = params.NB_VAR_ESTIMATORS !== 0 ? '--save-ic' : ''
        outprefix = "tmle." + parameterfile.getName().replace(".yaml", "")
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --threads=${task.cpus} --startup-file=no /TargetedEstimation.jl/scripts/tmle.jl \
        $data $parameterfile $estimatorfile $outprefix \
        $save_ic \
        --pval-threshold=${params.PVAL_SIEVE}
        """
}

process TMLEInputsFromParamFiles {
    container "olivierlabayle/tl-core:v0.3.0"
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
    container "olivierlabayle/tl-core:v0.3.0"
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
        batch_size = params.PHENOTYPES_BATCH_SIZE == 0 ? "" :  "--phenotype-batch-size ${params.PHENOTYPES_BATCH_SIZE}"
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        extra_covariates = extra_covariates.name != 'NO_EXTRA_COVARIATE' ? "--extra-covariates $extra_covariates" : ''
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        from-actors $bqtls $trans_actors_prefix $extra_confounders $extra_treatments $extra_covariates --orders ${params.ORDERS}
        """
}
