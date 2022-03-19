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

@testset "Test tmle_phenotypes_batches" begin
    n_phenotypes = 10
    temp_file = "temp_phenotypes.csv"
    phenotypes_names = ["x$i" for i in 1:n_phenotypes]
    CSV.write(
        temp_file, 
        DataFrame(rand(10, n_phenotypes+1), vcat("SAMPLE_ID", phenotypes_names))
        )

    # Test basic mode
    parsed_args = Dict(
        "phenotypes-file" => temp_file,
        "batch-size" => "3"
    )
    UKBBEpistasisPipeline.tmle_phenotypes_batches(parsed_args)

    expected_batches = [
        ["x1", "x2", "x3"], 
        ["x4", "x5", "x6"], 
        ["x7", "x8", "x9"], 
        ["x10"]
    ]
    for batch in 1:4
        batch_file = UKBBEpistasisPipeline.batch_filename(batch)
        @test readlines(open(batch_file)) == expected_batches[batch]
        rm(batch_file)
    end

    # Test max mode
    parsed_args = Dict(
        "phenotypes-file" => temp_file,
        "batch-size" => "max"
    )
    UKBBEpistasisPipeline.tmle_phenotypes_batches(parsed_args)

    batch_file = UKBBEpistasisPipeline.batch_filename(1)
    @test readlines(open(batch_file)) == phenotypes_names
    rm(batch_file)
    rm(temp_file)

end

end;

true