module TestGRM
using Test
using UKBBEpistasisPipeline

@testset "Test prepend_size" begin
    parsed_args = Dict(
        "in-grmfile"=> "data/grm/test.part_1000_0001.grm.bin", 
        "out-grmfile"=>"GRM_1.bin"
    )

    prepend_size(parsed_args)

    before_prepend = UKBBEpistasisPipeline.grm_part_from_gcta(parsed_args["in-grmfile"])
    after_prepend = UKBBEpistasisPipeline.grm_part(parsed_args["out-grmfile"])

    @test before_prepend == after_prepend

    rm(parsed_args["out-grmfile"])
end


end;

true