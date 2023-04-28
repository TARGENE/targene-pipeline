using Test
using CSV
using DataFrames
using YAML
using TMLE

# Move to the project's root
cd(dirname(dirname(@__FILE__)))

const SUMMARY_COLUMNS = [
    "PARAMETER_TYPE", "TREATMENTS", "CASE", "CONTROL", "TARGET", "CONFOUNDERS", "COVARIATES", 
    "INITIAL_ESTIMATE", "TMLE_ESTIMATE", "TMLE_STD", "TMLE_PVALUE", "TMLE_LWB", "TMLE_UPB", 
    "ONESTEP_ESTIMATE", "ONESTEP_STD", "ONESTEP_PVALUE", "ONESTEP_LWB", "ONESTEP_UPB", "LOG"
]

const SIEVE_COLUMNS = ["SIEVE_STD", "SIEVE_PVALUE", "SIEVE_LWB", "SIEVE_UPB"]

"""
    check_fails_are_extremely_rare_traits(output, dataset)

The failures are due to traits with only very few cases (<=3) which is due to the fact that 
the dataset is quite degenerate.
"""
function check_fails_are_extremely_rare_traits(output, dataset)
    groups = groupby(output, :TREATMENTS)
    for (treatment, group) in pairs(groups)
        treatment = treatment.TREATMENTS
        fails = filter(x -> x.LOG !== missing, group)
        for target in unique(fails.TARGET)
            @test eltype(dataset[!, target]) == Bool
            nomissing = dropmissing(dataset, Symbol.([target, split(treatment, "_&_")...]))
            @test sum(nomissing[!, target]) <= 3
        end
    end
end

"""

This is more of a non-regression test where we check that at least a certain amount of
estimates are a success.
"""
function test_n_success_more_than_threshold(output, threshold)
    successes = filter(x -> x.LOG === missing, output)
    n_succeses_per_treatment = DataFrames.combine(groupby(successes, :TREATMENTS),  nrow)
    @test all(x > threshold for x in n_succeses_per_treatment.nrow)
end
