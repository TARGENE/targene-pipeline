module TestNullSimulation

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test Null Simulation Workflow" begin
    cmd = `nextflow run main.nf -entry NULL_SIMULATION -c test/configs/simulation.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
end

end