module TestGRM
using Test
using UKBBEpistasisPipeline

@testset "Test readGRM" begin
    prefix = joinpath("data", "grm", "test.grm")
    GRM, ids = UKBBEpistasisPipeline.readGRM(prefix)
    @test size(GRM, 1) == 18915
    @test size(ids, 1) == 194
end

end;

true