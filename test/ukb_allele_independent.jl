module TestUKBAlleleIndependent

using Test
using JLD2
using Arrow
using DataFrames
using TMLE
using CSV
using TargetedEstimation
using Serialization

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test ukb_from_actors.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_allele_independent.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    ## Checking main output
    # Results
    hdf5_results_file = jldopen(joinpath("results", "results.hdf5"))
    results_from_hdf5 = vcat((hdf5_results_file[key] for key in keys(hdf5_results_file))...)
    results_from_json = TMLE.read_json(joinpath("results", "results.json"))
    @test length(results_from_json) == length(results_from_hdf5) > 300

    failed_results = retrieve_failed_results(results_from_hdf5; expected_keys=(:TMLE, :OSE))
    # All fails are due to fluctuation failure due non positive definite matrix
    # This does not affect the OSE
    @test isempty(failed_results.OSE)
    # Less than 1/3 affected: this is still quite significant
    @test length(failed_results.TMLE) / length(results_from_hdf5) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE)

    dataset = Arrow.Table(joinpath("results", "dataset.arrow")) |> DataFrame

    check_fails_are_extremely_rare_traits(failed_results.TMLE, dataset; ncases=3)

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    bQTLs = Symbol.(CSV.read(joinpath("test", "data", "actors", "bqtls.csv"), DataFrame).ID)

    config_1 = deserialize(joinpath("results", "estimands", "final.estimands_1.jls"))
    @test length(config_1.estimands) == 400
    config_2 = deserialize(joinpath("results", "estimands", "final.estimands_2.jls"))
    @test 0 < length(config_2.estimands) <= 400
    for Ψ in vcat(config_1.estimands, config_2.estimands)
        @test length(intersect(keys(Ψ.treatment_values), bQTLs)) == 1
    end
end

end