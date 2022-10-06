process UKBFieldsList {
    container "olivierlabayle/ukbmain:extract_subset_split"

    input:
        path traits_config

    output:
        path "fields_list.txt"

    script:
        traits_config = traits_config.name != 'NO_UKB_TRAIT_CONFIG' ? "--conf $traits_config" : ''
        "julia --project=/UKBMain.jl --startup-file=no /UKBMain.jl/scripts/build_fields_list.jl $traits_config"
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
        path 'traits.csv'

    script:
        traits_config = traits_config.name != 'NO_UKB_TRAIT_CONFIG' ? "--conf $traits_config" : ''
        withdrawal_list = withdrawal_list.name != 'NO_WITHDRAWAL_LIST' ? "--withdrawal-list $withdrawal_list" : ''
        """
        julia --project=/UKBMain.jl --startup-file=no /UKBMain.jl/scripts/process_main_dataset.jl \
        $dataset $traits_config $withdrawal_list
        """
}