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
    container "olivierlabayle/targeted-estimation:0.5"
    publishDir "$params.OUTDIR/csvs",  mode: 'symlink', pattern: "*.csv"
    publishDir "$params.OUTDIR/hdf5files/inf_curves",  mode: 'symlink', pattern: "*.hdf5"
    label "bigmem"
    label "multithreaded"

    input:
        path data
        path parameterfile
        path estimatorfile
    
    output:
        path "${csvout}", emit: tmle_csv
        path "${hdf5out}", optional: true, emit: inf_curve
    
    script:
        basename = "tmle." + parameterfile.getName().take(parameterfile.getName().lastIndexOf('.'))
        csvout = basename + ".csv"
        hdf5out = basename + ".hdf5"
        hdf5option = params.SAVE_IC == true ? "--hdf5-out=${hdf5out}" : ""
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --procs=${task.cpus} --threads=${task.cpus} --startup-file=no /TargetedEstimation.jl/scripts/tmle.jl \
        $data $parameterfile $estimatorfile $csvout \
        $hdf5option \
        --chunksize=100 \
        --pval-threshold=${params.PVAL_THRESHOLD}
        """
}

process TMLEInputsFromParamFile {
    container "olivierlabayle/tl-core:0.5"
    publishDir "$params.OUTDIR/parameters", mode: 'symlink', pattern: "*.yaml"
    publishDir "$params.OUTDIR/tmle_inputs", mode: 'symlink', pattern: "*.csv"
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path parameter

    output:
        path "final.data.csv", emit: traits
        path "final.*.yaml", emit: parameters

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        $batch_size \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        from-param-file $parameter
        """
}

process TMLEInputsFromActors {
    container "olivierlabayle/tl-core:0.5"
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
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        extra_covariates = extra_covariates.name != 'NO_EXTRA_COVARIATE' ? "--extra-covariates $extra_covariates" : ''
        genotypes_as_int = params.GENOTYPES_AS_INT == false ? "" : "--genotypes-as-int"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        from-actors $bqtls $trans_actors_prefix $extra_confounders $extra_treatments $extra_covariates --orders ${params.ORDERS} ${genotypes_as_int}
        """
}
