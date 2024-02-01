include { TMLEInputsFromActors; TMLEInputsFromParamFile } from '../modules/estimation_inputs.nf'

workflow EstimationInputs {
    take:
        study_design
        bgen_files
        traits
        genetic_confounders
        estimands_file
        bqtls_file
        transactors_files
        extra_confounders
        extra_treatments
        extra_covariates
        batch_size
        call_threshold
        positivity_constraint

    main:
        if (study_design == "FROM_ACTORS") {
            tmle_inputs = TMLEInputsFromActors(
                bgen_files,
                traits,
                genetic_confounders,
                extra_confounders,
                extra_treatments,
                extra_covariates,
                bqtls_file,
                transactors_files,
                batch_size,
                call_threshold,
                positivity_constraint
                )
        }
        else if (study_design == "CUSTOM"){
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                estimands_file,
                batch_size,
                call_threshold,
                positivity_constraint,
                "from-param-file"
                )
        }
        else if (study_design == "ALLELE_INDEPENDENT"){
            tmle_inputs = TMLEInputsFromParamFile(
                bgen_files,
                traits,
                genetic_confounders,
                estimands_file,
                batch_size,
                call_threshold,
                positivity_constraint,
                "allele-independent"
                )
        }
        else { 
            throw new Exception("This STUDY_DESIGN is not available.")
        }
    
    emit:
        aggregated_dataset = tmle_inputs.dataset
        estimands = tmle_inputs.estimands
}