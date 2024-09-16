module TestUKBAlleleIndependent

using Test
using JLD2
using TMLE
using TMLECLI
using Serialization
using YAML

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"] 

include("utils.jl")

@testset "Test ukb_interactions_group.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_interactions_group.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    # QQ plot
    @test isfile("results/QQ.png")

    # Check HDF5 results
    results = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    nresults = size(results, 1)
    @test nresults > 20
    nfails = 0
    treatment_combinations = Set([])
    for Ψ̂ in results.TMLE_GLM_GLM
        if Ψ̂ isa TMLECLI.FailedEstimate
            nfails += 1
        else
            Ψ = first(Ψ̂.estimand.args)
            push!(treatment_combinations, Tuple(keys(Ψ.treatment_values)))
        end
    end
    @test nfails / nresults < 1/4
    @test treatment_combinations ==  Set([ 
        (Symbol("2:14983:G:A"), Symbol("3:3502414:T:C")),
        (Symbol("1:238411180:T:C"), Symbol("2:14983:G:A"))
    ])

    # Check summary file
    summary_results = YAML.load_file(joinpath("results", "results.summary.yaml"))
    @test size(summary_results, 1) == nresults

    # Dataset
    dataset = TMLECLI.instantiate_dataset(joinpath("results", "datasets", "all_genotypes.data.arrow"))
    @test Set(names(dataset)) == Set(vcat("SAMPLE_ID", TRAITS, PCS, ["2:14983:G:A", "3:3502414:T:C", "1:238411180:T:C"]))

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end