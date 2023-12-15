include { longest_prefix } from './utils.nf'

process TMLEInputsFromParamFile {
    container "olivierlabayle/tl-core:cvtmle"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path genetic_confounders
        path parameter

    output:
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        batch_size = params.BATCH_SIZE == 0 ? "" : "--batch-size ${params.BATCH_SIZE}"
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
    container "olivierlabayle/tl-core:cvtmle"
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
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
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        trans_actors_prefix = longest_prefix(trans_actors)
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batch-size ${params.BATCH_SIZE}"
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        extra_covariates = extra_covariates.name != 'NO_EXTRA_COVARIATE' ? "--extra-covariates $extra_covariates" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/tmle_inputs.jl \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${params.CALL_THRESHOLD} \
        --pcs $genetic_confounders \
        --positivity-constraint ${params.POSITIVITY_CONSTRAINT} \
        $batch_size \
        from-actors $bqtls $trans_actors_prefix $extra_confounders $extra_treatments $extra_covariates --orders ${params.ORDERS}
        """
}

workflow EstimationInputs {
    take:
        traits
        genetic_confounders

    main:
        bgen_files = Channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect()

        if (params.STUDY_DESIGN == "FROM_ACTORS") {
            bqtls = Channel.value(file("$params.BQTLS"))
            trans_actors = Channel.fromPath("$params.TRANS_ACTORS", checkIfExists: true).collect()
            extra_confounders = Channel.value(file("$params.EXTRA_CONFOUNDERS"))
            extra_treatments = Channel.value(file("$params.ENVIRONMENTALS"))
            extra_covariates = Channel.value(file("$params.EXTRA_COVARIATES"))
            tmle_inputs = TMLEInputsFromActors(
                bgen_files,
                traits,
                genetic_confounders,
                extra_confounders,
                extra_treatments,
                extra_covariates,
                bqtls,
                trans_actors)
        }
        else if (params.STUDY_DESIGN == "CUSTOM"){
            estimands_file = Channel.value(file("$params.ESTIMANDS_FILE"))
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                estimands_file)
        }
        else { 
            throw new Exception("This STUDY_DESIGN is not available.")
        }
    
    emit:
        aggregated_dataset = tmle_inputs.dataset
        estimands = tmle_inputs.estimands
}