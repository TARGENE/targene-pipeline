module TestGRM
using Test
using UKBBEpistasisPipeline
using CSV
using DataFrames
using Arrow

@testset "Test GRMfromGCTAFile" begin
    parsed_args = Dict("inprefix"=> "data/grm/test", "outprefix"=>"GRM")
    arrow_file = "GRM.arrow"
    ids_file = "GRM.ids.csv"
    grm_parts_to_arrow(parsed_args)

    arrow_grm = Arrow.Table(arrow_file).RELATIONSHIPS
    ids = CSV.File(ids_file) |> DataFrame

    # This "original" file has been generated using 
    # the provided R code from the gcta website
    expected_grm = CSV.File(parsed_args["inprefix"]*".grm.txt", header=false) |> DataFrame
    expected_ids = CSV.File(parsed_args["inprefix"]*".grm.id", header=[:FAMILY_ID, :ID]) |> DataFrame
    
    @test arrow_grm ≈ expected_grm[!, 1]
    @test ids == expected_ids
    rm(arrow_file)
    rm(ids_file)
end


end;

true