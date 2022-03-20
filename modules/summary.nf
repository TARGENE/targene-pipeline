process Summary {
    input:
        tuple val(rsids), file(tmle), file(sieve)
    
    output:
        path "${rsids}_summary.csv"
    
    script:
        "touch ${rsids}_summary.csv"
}