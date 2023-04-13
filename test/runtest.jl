using Test
using CSV
using DataFrames
using StatsBase
using YAML

# Move to the project's root
cd(dirname(pwd()))

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


@testset "Test from_actors.config" begin
    r = run(`nextflow run main.nf -c conf/ci_jobs/from_actors.config -profile eddie -resume`)
    @test r.exitcode == 0

    ## Checking main output
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = CSV.read(joinpath("results", "tmle_inputs", "final.data.csv"), DataFrame)

    @test names(output) == vcat(SUMMARY_COLUMNS, SIEVE_COLUMNS)
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["1:238411180:T:C", "3:3502414:T:C", "1:238411180:T:C_&_2:14983:G:A", "3:3502414:T:C_&_2:14983:G:A"])
    
    check_fails_are_extremely_rare_traits(output, dataset)
    test_n_success_more_than_threshold(output, 20)

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_1.yaml"))
    @test param_dict["T"] == ["1:238411180:T:C"]
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_2.yaml"))
    @test param_dict["T"] == ["3:3502414:T:C"]
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_3.yaml"))
    @test param_dict["T"] == ["1:238411180:T:C", "2:14983:G:A"]
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_4.yaml"))
    @test param_dict["T"] == ["3:3502414:T:C", "2:14983:G:A"]
end


@testset "Test from_actors.config" begin
    r = run(`nextflow run main.nf -c conf/ci_jobs/from_param_files.config -profile eddie -resume`)
    @test r.exitcode == 0
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = CSV.read(joinpath("results", "tmle_inputs", "final.data.csv"), DataFrame)
    @test names(output) == SUMMARY_COLUMNS
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["3:3502414:T:C_&_1:238411180:T:C", "2:14983:G:A"])
    
    check_fails_are_extremely_rare_traits(output, dataset)
    test_n_success_more_than_threshold(output, 20)
end