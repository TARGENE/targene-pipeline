default_args = ["eddie", "-resume"] 
for (i, arg) in enumerate(ARGS)
    default_args[i] = arg
end

include("utils.jl")

@testset "Test from_param_files.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/from_param_files.config $default_args`
    @info string("The following command will be run:\n", cmd)
    @test r.exitcode == 0
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = CSV.read(joinpath("results", "tmle_inputs", "final.data.csv"), DataFrame)
    @test names(output) == SUMMARY_COLUMNS
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["3:3502414:T:C_&_1:238411180:T:C", "2:14983:G:A"])
    
    check_fails_are_extremely_rare_traits(output, dataset)
    test_n_success_more_than_threshold(output, 20)
end
