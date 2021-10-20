module TestPhenotypes

using Test
using UKBBEpistasisPipeline
using CSV
using DataFrames


@testset "Test prepare_phenotypes" begin
    # Test with no phenotype list
    parsed_args = Dict(
        "binary-phenotypes" => joinpath("data", "ukbb", "binary_phenotypes.csv"),
        "continuous-phenotypes" => joinpath("data", "ukbb", "quantitative_phenotypes.csv"),
        "bridge" => joinpath("data", "ukbb", "bridge.csv"),
        "withdrawal-list" => joinpath("data", "ukbb", "withdrawal_list.csv"),
        "output" => "full_phenotypes.csv",
        "phenotypes-list" => nothing
    )
    UKBBEpistasisPipeline.prepare_phenotypes(parsed_args)
    
    phenotypes = CSV.File(parsed_args["output"]) |> DataFrame

    @test names(phenotypes) == ["eid",
                                "geneatlas_id",
                                "clinical_c_Block_J40-J47",
                                "clinical_c_G25",
                                "clinical_c_G20",
                                "1408-0.0",
                                "1777-0.0",
                                "1727-0.0",
                                "1548-0.0"]
    @test phenotypes.eid == [1000023, 1000030]
    @test phenotypes.geneatlas_id == [1048581, 5242887]

    rm(parsed_args["output"])

    # Test with phenotypes list, no withdrawal-list

    parsed_args = Dict(
        "binary-phenotypes" => joinpath("data", "ukbb", "binary_phenotypes.csv"),
        "continuous-phenotypes" => joinpath("data", "ukbb", "quantitative_phenotypes.csv"),
        "bridge" => joinpath("data", "ukbb", "bridge.csv"),
        "output" => "full_phenotypes.csv",
        "withdrawal-list" => nothing,
        "phenotypes-list" => "clinical_c_G25,clinical_c_G20,1408-0.0"
    )
    UKBBEpistasisPipeline.prepare_phenotypes(parsed_args)
    
    phenotypes = CSV.File(parsed_args["output"]) |> DataFrame

    @test names(phenotypes) == ["eid",
                                "geneatlas_id",
                                "clinical_c_G25",
                                "clinical_c_G20",
                                "1408-0.0"]
    @test phenotypes.eid == [1000014, 1000023, 1000030]
    @test phenotypes.geneatlas_id == [2097158, 1048581, 5242887]

    rm(parsed_args["output"])

end

end;

true