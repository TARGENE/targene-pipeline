
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

process FromASBxTransActors {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/queries", mode: 'symlink'
    label "bigmem"

    input:
        path bgenfiles
        path asbs
        path trans_actors
        path excluded_snps
        path sample_ids

    output:
        path "outputs/*.toml", emit: queries
        path "outputs/genotypes.csv", emit: genotypes

    script:
        chr_prefix = longest_prefix(bgenfiles)
        asb_prefix = longest_prefix(asbs)
        exclude = excluded_snps.name != 'NO_FILE' ? "--exclude $excluded_snps" : ''
        """
        mkdir -p outputs
        julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/build_genotypes_and_queries.jl \
        $chr_prefix --mode=asb --outdir=outputs --call-threshold=${params.CALL_THRESHOLD} $exclude \
        --minor-genotype-freq=${params.MINOR_CAT_FREQUENCY} --asb-prefix=$asb_prefix --trans-actors=$trans_actors --sample-ids=$sample_ids
        """
}

process FromGivenQueries {
    container "olivierlabayle/tl-core:sample_filtering"
    publishDir "$params.OUTDIR/queries", mode: 'symlink'
    label "bigmem"

    input:
        path bgenfiles
        path query_files
        path excluded_snps
        path sample_ids

    output:
        path "outputs/*.toml", emit: queries
        path "outputs/genotypes.csv", emit: genotypes

    script:
        chr_prefix = longest_prefix(bgenfiles)
        query_prefix = longest_prefix(query_files)
        exclude = excluded_snps.name != 'NO_FILE' ? "--exclude $excluded_snps" : ''
        """
        mkdir -p outputs
        julia --project=/TMLEEpistasis.jl --startup-file=no /TMLEEpistasis.jl/bin/build_genotypes_and_queries.jl \
        $chr_prefix --mode=given --outdir=outputs --call-threshold=${params.CALL_THRESHOLD} $exclude \
        --minor-genotype-freq=${params.MINOR_CAT_FREQUENCY} --query-prefix=$query_prefix --sample-ids=$sample_ids
        """
}
