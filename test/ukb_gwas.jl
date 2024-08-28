module TestFromParamFile

using Test
using DataFrames
using TMLE
using Arrow
using JLD2
using TMLECLI
using TargeneCore

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

@testset "Test ukb_estimands_file.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_gwas.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    # Check datasets are consistent
    traits_and_pcs = ["SAMPLE_ID", "Body mass index (BMI)", "Number of vehicles in household", "Cheese intake", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]
    for chr in 1:3
        dataset = DataFrame(Arrow.Table(joinpath("results", "datasets", "ukb_chr$(chr).data.arrow")))
        @test nrow(dataset) == 500
        columnnames = names(dataset)
        @test issubset(traits_and_pcs, columnnames)
        variant_columns = setdiff(columnnames, traits_and_pcs)
        @test all(startswith(colname, string(chr)) for colname in variant_columns)
    end

    # Check results
    results = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    @test length(results) > 80
    variants = []
    failed_results = []
    for estimators_results in results
        for (estimatoir, Ψ̂) in zip(keys(estimators_results), estimators_results)
            if Ψ̂ isa TMLECLI.FailedEstimate
                push!(failed_results, Ψ̂)
            end
            Ψ = Ψ̂.estimand
            @test Ψ isa JointEstimand
            @test get_outcome(Ψ) == Symbol("Body mass index (BMI)")
            variant = only(TargeneCore.get_treatments(Ψ))
            push!(variants, variant)
            @test get_confounders(Ψ, variant) == (:PC1, :PC2, :PC3, :PC4, :PC5, :PC6)
            @test get_outcome_extra_covariates(Ψ) == Symbol.(("Cheese intake", "Number of vehicles in household"))
        end
    end
    @test Set(string(v)[1] for v ∈ variants) == Set(['1', '2', '3'])
    @test length(failed_results) == 0

    # Check QQ
    @test isfile(joinpath("results", "QQ.png"))

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 100
end

end

true