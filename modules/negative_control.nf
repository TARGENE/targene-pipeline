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

process GeneratePermutationTestsData {
    container "olivierlabayle/negative-controls:cvtmle"
    publishDir "${params.OUTDIR}/permutation_tests", mode: 'symlink', pattern: '*.arrow'
    publishDir "${params.OUTDIR}/permutation_tests/estimands", mode: 'symlink', pattern: '*.jls'
    label "bigmem"
    
    input:
        path dataset
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.jls", emit: estimands

    script:
        limit = params.MAX_PERMUTATION_TESTS == null ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_permutation_data.jl \
        ${dataset} ${results} \
        ${limit} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --orders=${params.PERMUTATION_ORDERS} \
        --chunksize=${params.BATCH_SIZE} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

process GenerateRandomVariantsTestsData {
    container "olivierlabayle/negative-controls:cvtmle"
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"
    
    input:
        path trans_actors
        path bgenfiles
        path results

    output:
        path "random_variants_estimands.jls"

    script:
        trans_actors_prefix = longest_prefix(trans_actors)
        bgen_prefix = longest_prefix(bgenfiles)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_random_variant_parameters.jl \
        ${trans_actors_prefix} ${results} ${bgen_prefix} \
        --out=random_variants_estimands.jls \
        --p=${params.N_RANDOM_VARIANTS} \
        --reltol=${params.MAF_MATCHING_RELTOL} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}