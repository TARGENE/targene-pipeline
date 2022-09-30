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
        path 'traits.csv'

    script:
        """
        julia --project=/UKBMain.jl --startup-file=no /UKBMain.jl/scripts/process_main_dataset.jl \
        $dataset --conf $traits_config --withdrawal-list $withdrawal_list --out traits.csv
        """
}