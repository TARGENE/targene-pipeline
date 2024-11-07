module TestRealisticSimulation

using Test
using TMLECLI
using JLD2
using DataFrames

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test Realistic Simulation Workflow" begin
    cmd = `nextflow run main.nf -entry REALISTIC_SIMULATION -c test/configs/simulation.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    results = jldopen(io -> io["results"], joinpath("results", "realistic_simulation_results.hdf5"))

    @test all(length(x) == 4 - nf for (x, nf) in zip(results.ESTIMATES, results.N_FAILED)) # 2 bootstraps per run * 2 random seeds
    @test all(x == 1000 for x in results.SAMPLE_SIZE)
    @test all(x -> x in [:OSE_GLM_GLM, :TMLE_GLM_GLM], results.ESTIMATOR)
    # For the first joint estimand: 
    # Total number of traits = 11 - (Number of vehicles in household + Skin colour) = 9
    # For the second estimand, only 1 outcome
    # Total = 9 + 1 = 10 estimands per-estimator 
    # 2 estimators so 20 total estimands
    @test nrow(results) == 20

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000

end

end