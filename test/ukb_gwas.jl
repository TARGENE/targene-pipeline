module TestFromParamFile

using Test
using DataFrames
using TMLE
using Arrow
using JLD2
using TmleCLI

# "local" profile assumes singularity is installed
args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

include("utils.jl")

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
    variants = []
    for estimators_results in results
        for (estimatoir, result) in zip(keys(estimators_results), estimators_results)
            push!(variants, TmleCLI.get_treatments(result.estimand))
        end
    end
    @test length(results) > 40

    failed_results = retrieve_failed_results(results)

    # All fails are due to fluctuation failure due non positive definite matrix
    # This does not affect the OSE
    @test isempty(failed_results.OSE_GLM_GLM)
    # Less than 1/3 affected: this is still quite significant
    @test length(failed_results.TMLE_GLM_GLM) / length(results) < 1/3
    @test all(startswith(x.msg, "Could not fluctuate") for x ∈ failed_results.TMLE_GLM_GLM)

    check_fails_are_extremely_rare_traits(failed_results.TMLE_GLM_GLM, dataset; ncases=3)

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 1000
end

end

true