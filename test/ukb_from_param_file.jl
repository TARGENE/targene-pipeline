module TestFromParamFile

using Test
using CSV
using DataFrames
using YAML
using TMLE
using Arrow
using Serialization
using JLD2

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

include(joinpath(@__DIR__, "test", "utils.jl"))

@testset "Test ukb_from_param_files.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/ukb_from_param_file.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    result_file = jldopen(joinpath("results", "results.hdf5"))
    results = vcat(results["Batch_1"], result_file["Batch_2"])
    dataset = DataFrame(Arrow.Table(joinpath("results", "dataset.arrow")))

    failed_results = (TMLE = [], OSE = [])
    for result ∈ results
        @test keys(result) == (:TMLE, :OSE)
        @test result.TMLE isa Union{TMLE.TMLEstimate, TargetedEstimation.FailedEstimate}
        @test result.OSE isa Union{TMLE.OSEstimate, TargetedEstimation.FailedEstimate}
        if result.TMLE isa TargetedEstimation.FailedEstimate
            push!(failed_results.TMLE, result.TMLE)
        end
        if result.OSE isa TargetedEstimation.FailedEstimate
            push!(failed_results.OSE, result.OSE)
        end
    end
    # All fails are due to fluctuation failure due non positive definite matrix
    # This does not affect the OSE
    @test isempty(failed_results.OSE)
    # Less than 1/3 affected: this is still quite significant
    @test length(failed_results.TMLE) / length(results) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE)

    check_fails_are_extremely_rare_traits(failed_results.TMLE, dataset; ncases=3)
end

end

true