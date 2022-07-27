process UKBFieldsList {
    container "olivierlabayle/ukbmain:extract_subset_split"

    input:
        path traits_config

    output:
        path "fields_list.txt"

    script:
        "julia --project=/UKBMain.jl --startup-file=no /UKBMain.jl/scripts/build_fields_list.jl --conf $traits_config"
}

process UKBConv {
    container "olivierlabayle/ukbmain:extract_subset_split"

    input:
        path fields_list
        path encrypted_dataset
        path encoding_file

    output:
        path "ukb_dataset.csv"

    script:
        "ukbconv $encrypted_dataset csv -i$fields_list -e$encoding_file -oukb_dataset"
}

process TraitsFromUKB {
    publishDir "$params.OUTDIR/traits", mode: 'symlink'
    container "olivierlabayle/ukbmain:extract_subset_split"
    memory "20G"

    input:
        path dataset
        path traits_config
        path withdrawal_list
    
    output:
        path 'processed.sample_ids.txt', emit: sample_ids
        path 'processed.continuous.phenotypes.csv', emit: continuous_phenotypes, optional: true
        path 'processed.binary.phenotypes.csv', emit: binary_phenotypes, optional: true
        path 'processed.covariates.csv', emit: covariates, optional: true
        path 'processed.extra_confounders.csv', emit: confounders, optional: true
        path 'processed.extra_treatments.csv', emit: treatments, optional: true
    
    script:
        """
        julia --project=/UKBMain.jl --startup-file=no /UKBMain.jl/scripts/process_main_dataset.jl \
        $dataset --conf $traits_config --withdrawal-list $withdrawal_list --out-prefix processed
        """
}