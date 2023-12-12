module TestFromParamFile

using Test
using CSV
using DataFrames
using YAML
using TMLE
using Arrow
using Serialization

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

include(joinpath(@__DIR__, "test", "utils.jl"))

@testset "Test ukb_from_param_files.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/ukb_from_param_file.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = DataFrame(Arrow.Table(joinpath("results", "tmle_inputs", "final.data.arrow")))
    @test names(output) == vcat(SUMMARY_COLUMNS, ADJUTMENT_COL)
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["1:238411180:T:C_&_3:3502414:T:C", "2:14983:G:A"])
    
    check_fails_are_extremely_rare_traits(output, dataset)
    test_n_success_more_than_threshold(output, 20)
end

end

true