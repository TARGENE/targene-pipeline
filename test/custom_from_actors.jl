# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test custom_from_actors.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/custom_from_actors.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = DataFrame(Arrow.Table(joinpath("results", "tmle_inputs", "final.data.arrow")))
    bQTLs = Symbol.(CSV.read(joinpath("test", "data", "actors", "bqtls.csv"), DataFrame).ID)

    @test names(output) == vcat(SUMMARY_COLUMNS, SIEVE_COLUMNS, ADJUTMENT_COL)
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["1:238411180:T:C", "3:3502414:T:C", "1:238411180:T:C_&_2:14983:G:A", "3:3502414:T:C_&_2:14983:G:A"])
        
    test_n_success_more_than_threshold(output, 20)

    # Here we test that the process generateIIDGenotypes has been run once for each chromosome
    n_chr_files = filter(x -> startswith(x, "LDpruned.filtered.ukb_chr"), readdir(joinpath("results","ld_pruned_chromosomes")))
    # There should be 4 files for each chromosome
    @test length(n_chr_files) == 4*3

    ## Checking parameter files correspond to either bQTL only or bQTL/eQTL
    parameters_tf1 = parameters_from_yaml(joinpath("results", "parameters", "final.TF1.param_1.yaml"))
    parameters_tf2 = parameters_from_yaml(joinpath("results", "parameters", "final.TF2.param_1.yaml"))
    @test size(parameters_tf1, 1) == 296
    @test size(parameters_tf2, 1) == 148
    for Ψ in vcat(parameters_tf1, parameters_tf2)
        @test keys(Ψ.treatment)[1] ∈ bQTLs
    end

end
