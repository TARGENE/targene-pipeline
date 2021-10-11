module TestFilterASB

using Test
using CSV
using DataFrames
using UKBBEpistasisPipeline

@testset "Test default_filter_asb" begin
    asb_files = ["data/asb_files/asb_1.csv", "data/asb_files/asb_2.csv"]
    out = joinpath("data", "asb_files", "asb-filtered-test.csv")
    # Test default mode
    parsed_args = Dict(
        "mode" => "default",
        "out" => out,
        "asb-files" => asb_files)
    filter_asb(parsed_args)

    filtered_asb = CSV.File(out) |> DataFrame
    @test filtered_asb.ID == ["RSID_2", "RSID_102"]

    rm(out)

    # Test unknown mode throws
    parsed_args = Dict(
        "mode" => "custom",
        "out" => out,
        "asb-files" => asb_files)
    @test_throws ArgumentError filter_asb(parsed_args)
    
end

end;

true