module TestUKBAlleleIndependent

using Test
using Arrow
using DataFrames

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test MAKE_DATASET" begin
    cmd = `nextflow run main.nf -entry MAKE_DATASET -c test/configs/dataset.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    dataset = Arrow.Table("results/dataset.arrow") |> DataFrame

end

end

true