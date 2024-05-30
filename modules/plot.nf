process GenerateSummaryPlots {
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    label 'targenecore_image'

    input:
        path resultsfile

    output:
        path "*.png"

    script:
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/bin/generate_summary_plots.jl \
        ${resultsfile} \
        --verbosity=${params.VERBOSITY} \
        --out-prefix="."
        """
}