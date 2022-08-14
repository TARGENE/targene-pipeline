process Summary {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/summaries", mode: 'symlink'
    label "bigmem"

    input:
        tuple val(treatment_id), file(tmle), file(sieve)
    
    output:
        path "${treatment_id}.summary.csv"
    
    script:
        sieve_option = sieve.getName() != 'NO_FILE' ? "--sieve" : ''
        """
        julia --project=/TargeneCore.jl --startup-file=no /TargeneCore.jl/bin/summarize.jl \
        final ${treatment_id}.summary.csv $sieve_option
        """
}