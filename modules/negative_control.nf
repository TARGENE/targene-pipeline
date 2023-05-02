process GeneratePermutationTestsData {
    container "olivierlabayle/negative-controls:initial_pipeline"
    publishDir "${params.OUTDIR}/permutation_data", mode: 'symlink'
    label "bigmem"
    
    input:
        path dataset
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.bin", emit: parameters

    script:
        limit = params.MAX_PERMUTATION_TESTS == null ? "" : "--limit=${params.MAX_PERMUTATION_TESTS}"
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_permutation_data.jl \
        ${dataset} ${results} \
        ${limit} \
        --pval-col=${params.PVAL_COL} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --orders=${params.PERMUTATION_ORDERS} \
        --chunksize=${params.BATCH_SIZE} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}

process GenerateRandomVariantsTestsData {
    container "olivierlabayle/negative-controls:initial_pipeline"
    publishDir "${params.OUTDIR}/random_variants_parameters", mode: 'symlink'
    label "bigmem"
    
    input:
        path trans_actors
        path bgenfiles
        path results

    output:
        path "permutation_dataset.arrow", emit: dataset
        path "*.bin", emit: parameters

    script:
        trans_actors_prefix = longest_prefix(trans_actors)
        bgen_prefix = longest_prefix(bgenfiles)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/NegativeControl  --startup-file=no /NegativeControl/bin/generate_random_variant_parameters.jl \
        ${trans_actors_prefix} ${results} ${bgen_prefix} \
        --p=${params.N_RANDOM_VARIANTS} \
        --reltol=${params.MAF_MATCHING_RELTOL} \
        --pval-col=${params.PVAL_COL} \
        --pval-threshold=${params.PVAL_THRESHOLD} \
        --chunksize=${params.BATCH_SIZE} \
        --rng=${params.RNG} \
        --verbosity=${params.VERBOSITY}
        """
}