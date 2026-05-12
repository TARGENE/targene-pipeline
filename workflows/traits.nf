include { ExtractTraits } from '../subworkflows/extract_traits.nf'

workflow EXTRACT_TRAITS {
    main:
        ukb_encoding_file = params.UKB_ENCODING_FILE
        ukb_config = channel.value(file("$params.UKB_CONFIG", checkIfExists: true))
        ukb_withdrawal_list = channel.value(file("$params.UKB_WITHDRAWAL_LIST", checkIfExists: true))
        traits_dataset = channel.value(file("$params.TRAITS_DATASET", checkIfExists: true))

        ExtractTraits(
            traits_dataset,
            ukb_config,
            ukb_withdrawal_list,
            ukb_encoding_file,
        )

    emit:
        ExtractTraits.out
}