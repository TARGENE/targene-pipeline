module TestCustomCohortFlatConfig

using Test
using Arrow
using DataFrames

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

@testset "Test custom_cohort_flat.config" begin
    cmd = `nextflow run main.nf -c test/configs/custom_cohort_flat.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end

true