include { longest_prefix } from './utils.nf'

process EstimationInputs {
    publishDir "$params.OUTDIR/estimands", mode: 'symlink', pattern: "*.jls"
    publishDir "$params.OUTDIR", mode: 'symlink', pattern: "*.arrow", saveAs: { filename -> "${params.ARROW_OUTPUT}" }
    label "bigmem"
    label 'targenecore_image'

    input:
        tuple val(genotypes_id), path(pcs), path(genotypes)
        path traits
        path config_file

    output:
        path "final.data.arrow", emit: dataset
        path "final.*.jls", emit: estimands

    script:
        genotypes_prefix = longest_prefix(genotypes)
        batch_size = params.BATCH_SIZE == 0 ? "" :  "--batchsize ${params.BATCH_SIZE}"
        call_threshold = params.CALL_THRESHOLD == null ? "" : "--call-threshold ${params.CALL_THRESHOLD}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl \
        estimation-inputs ${config_file} \
        --genotypes-prefix=${genotypes_prefix} \
        --traits-file=${traits} \
        --pcs-file=${pcs} \
        --outprefix=final \
        ${batch_size} \
        ${call_threshold} \
        --positivity-constraint=${params.POSITIVITY_CONSTRAINT} \
        --verbosity=${params.VERBOSITY}
        """
}