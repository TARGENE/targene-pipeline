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
        val batch_size
        val call_threshold
        val positivity_constraint
        val command

    output:
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        batchsize = batch_size == 0 ? "" : "--batch-size ${batch_size}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint ${positivity_constraint} \
        $batchsize \
        --out-prefix=final \
        --verbosity=${params.VERBOSITY} \
        $command $parameter \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${call_threshold} \
        --pcs $genetic_confounders \
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
        val batch_size
        val call_threshold
        val positivity_constraint

    output:
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        bgen_prefix = longest_prefix(bgenfiles)
        trans_actors_prefix = longest_prefix(trans_actors)
        batch_size = batch_size == 0 ? "" :  "--batch-size ${batch_size}"
        extra_confounders = extra_confounders.name != 'NO_EXTRA_CONFOUNDER' ? "--extra-confounders $extra_confounders" : ''
        extra_treatments = extra_treatments.name != 'NO_EXTRA_TREATMENT' ? "--extra-treatments $extra_treatments" : ''
        extra_covariates = extra_covariates.name != 'NO_EXTRA_COVARIATE' ? "--extra-covariates $extra_covariates" : ''
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_tl_inputs.jl \
        --positivity-constraint ${positivity_constraint} \
        $batch_size \
        --out-prefix=final \
        --verbosity=${params.VERBOSITY} \
        from-actors $bqtls $trans_actors_prefix \
        $extra_confounders \
        $extra_treatments \
        $extra_covariates \
        --orders ${params.ORDERS} \
        --traits $traits \
        --bgen-prefix $bgen_prefix \
        --call-threshold ${call_threshold} \
        --pcs $genetic_confounders \
        """
}