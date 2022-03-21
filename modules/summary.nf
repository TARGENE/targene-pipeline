process Summary {
    container "olivierlabayle/tmle-epistasis:0.3.0"
    publishDir "$params.OUTDIR/summaries", mode: 'symlink'
    label "bigmem"

    input:
        tuple val(rsids), file(tmle), file(sieve)
    
    output:
        path "${rsids}_summary.csv"
    
    script:
        """
        julia --project=/TMLEEpistasis.jl --startup-file=no --sysimage /TMLEEpistasis.jl/test/TMLEEpistasisSysimage.so \
        /TMLEEpistasis.jl/bin/summarize.jl $rsids ${rsids}_summary.csv
        """
}