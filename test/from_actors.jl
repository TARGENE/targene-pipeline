args = length(ARGS) > 0 ? ARGS : ["eddie", "-resume"] 

include("utils.jl")

@testset "Test from_actors.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/from_actors.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    ## Checking main output
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = CSV.read(joinpath("results", "tmle_inputs", "final.data.csv"), DataFrame)

    @test names(output) == vcat(SUMMARY_COLUMNS, SIEVE_COLUMNS)
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["1:238411180:T:C", "3:3502414:T:C", "1:238411180:T:C_&_2:14983:G:A", "3:3502414:T:C_&_2:14983:G:A"])
    
    check_fails_are_extremely_rare_traits(output, dataset)
    test_n_success_more_than_threshold(output, 20)

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_1.yaml"))
    @test param_dict["T"] == ["1:238411180:T:C"]
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_2.yaml"))
    @test param_dict["T"] == ["3:3502414:T:C"]
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_3.yaml"))
    @test param_dict["T"] == ["1:238411180:T:C", "2:14983:G:A"]
    param_dict = YAML.load_file(joinpath("results", "parameters", "final.param_4.yaml"))
    @test param_dict["T"] == ["3:3502414:T:C", "2:14983:G:A"]
end