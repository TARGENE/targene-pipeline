# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test custom_from_param_file.config" begin
    cmd = `nextflow run main.nf -c conf/ci_jobs/custom_from_param_file.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    output = CSV.read(joinpath("results", "summary.csv"), DataFrame)
    dataset = DataFrame(Arrow.Table(joinpath("results", "tmle_inputs", "final.data.arrow")))
    @test names(output) == vcat(SUMMARY_COLUMNS, ADJUTMENT_COL)
    # 2 bQTLs and 1 trans-actor
    @test Set(unique(output.TREATMENTS)) == Set(["1:238411180:T:C_&_3:3502414:T:C", "2:14983:G:A"])
    
    test_n_success_more_than_threshold(output, 20)

    # Here we test that the process generateIIDGenotypes has been run once for each chromosome
    n_chr_files = filter(x -> startswith(x, "LDpruned.filtered.ukb_chr"), readdir(joinpath("results","ld_pruned_chromosomes")))
    # There should be 4 files for each chromosome
    @test length(n_chr_files) == 4*3
end
