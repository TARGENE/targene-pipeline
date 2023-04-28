# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test from_actors.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/from_actors.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    ## Checking main output
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = CSV.read(joinpath("results", "tmle_inputs", "final.data.csv"), DataFrame)
    bQTLs = Symbol.(CSV.read(joinpath("test", "data", "actors", "bqtls.csv"), DataFrame).ID)

    @test names(output) == vcat(SUMMARY_COLUMNS, SIEVE_COLUMNS)
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["1:238411180:T:C", "3:3502414:T:C", "1:238411180:T:C_&_2:14983:G:A", "3:3502414:T:C_&_2:14983:G:A"])
    
    check_fails_are_extremely_rare_traits(output, dataset)
    test_n_success_more_than_threshold(output, 20)

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    parameters_1 = parameters_from_yaml(joinpath("results", "parameters", "final.param_1.yaml"))
    @test size(parameters_1, 1) == 400
    for Ψ in parameters_1
        @test keys(Ψ.treatment)[1] ∈ bQTLs
    end

    parameters_2 = parameters_from_yaml(joinpath("results", "parameters", "final.param_2.yaml"))
    @test size(parameters_2, 1) == 44
    for Ψ in parameters_2
        @test keys(Ψ.treatment)[1] ∈ bQTLs
    end

end