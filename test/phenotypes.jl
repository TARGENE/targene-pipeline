module TestPhenotypes

using Test
using UKBBEpistasisPipeline
using CSV
using DataFrames


@testset "Test prepare_phenotypes" begin
    # Test for continuous phenotypes with withdrawal list
    parsed_args = Dict(
        "phenotypes-file" => joinpath("data", "ukbb", "quantitative_phenotypes.csv"),
        "bridge" => joinpath("data", "ukbb", "bridge.csv"),
        "withdrawal-list" => joinpath("data", "ukbb", "withdrawal_list.csv"),
        "output" => "full_phenotypes.csv",
        "phenotypes-list" => nothing
    )
    UKBBEpistasisPipeline.prepare_phenotypes(parsed_args)
    
    phenotypes = CSV.File(parsed_args["output"]) |> DataFrame

    @test names(phenotypes) == ["SAMPLE_ID",
                                "1408-0.0",
                                "1777-0.0",
                                "1727-0.0",
                                "1548-0.0"]
    @test phenotypes.SAMPLE_ID == [1000023, 1000030]

    rm(parsed_args["output"])

    # Test for continuous phenotypes with withdrawal phenotypes list

    parsed_args = Dict(
        "phenotypes-file" => joinpath("data", "ukbb", "binary_phenotypes.csv"),
        "bridge" => joinpath("data", "ukbb", "bridge.csv"),
        "output" => "full_phenotypes.csv",
        "withdrawal-list" => nothing,
        "phenotypes-list" => joinpath("data", "ukbb", "phenotypes_list.csv")
    )
    UKBBEpistasisPipeline.prepare_phenotypes(parsed_args)
    
    phenotypes = CSV.File(parsed_args["output"]) |> DataFrame

    @test names(phenotypes) == ["SAMPLE_ID",
                                "clinical_c_G25",
                                "clinical_c_G20"]
    @test phenotypes.SAMPLE_ID == [1000190, 1000023, 1000014, 1000030]

    rm(parsed_args["output"])

end

end;

true