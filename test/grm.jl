module TestGRM
using Test
using UKBBEpistasisPipeline
using Mmap
using CSV
using DataFrames

@testset "Test grm_from_gcta" begin
    parsed_args = Dict(
        "inprefix"=> "data/grm/test", 
        "outprefix"=>"GRM"
    )

    @time grm_from_gcta(parsed_args)
    grm, ids = UKBBEpistasisPipeline.readGRM("GRM")
    
    expected_grm = CSV.File("data/grm/test.grm.txt", header = ["VALUE"]) |> DataFrame
    index = 1
    for i in 1:194
        for j in 1:i
            @test grm[i, j] == grm[j, i] ≈ expected_grm.VALUE[index]
            index += 1
        end
    end
    rm("GRM.bin")
    rm("GRM.ids.csv")
end


end;

true