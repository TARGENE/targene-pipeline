module TestCustomCohortFlatConfig

using Test
using Arrow
using DataFrames
using JLD2
using TMLECLI
using CSV
using TargeneCore
using YAML

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test that pipeline finishes even when no test can be computed" begin
    cmd = `nextflow run main.nf -c test/configs/empty_qq.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    dataset_origin = CSV.read(joinpath("test", "assets", "traits_sub.csv"), DataFrame)
    dataset = DataFrame(Arrow.Table(joinpath("results", "datasets", "all_genotypes.data.arrow")))
    @test issubset(names(dataset_origin), names(dataset))

    # Check QQ plot is still created
    @test isfile(joinpath("results", "QQ.png"))
    
    # Check full HDF5 results
    results = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    @test size(results, 1) == 1

    # Check summary file
    summary_results = YAML.load_file(joinpath("results", "results.summary.yaml"))
    @test size(summary_results, 1) == size(results, 1)
    
    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end

true