module TestSubsetMode

using Test
using DataFrames
using TMLE
using Arrow
using JLD2
using TMLECLI
using TargeneCore
using YAML

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

@testset "Test ukb_gwas_subset.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_gwas_subset.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0
    
    # Check datasets are created from subset BED files
    traits_and_pcs = ["SAMPLE_ID", "Body mass index (BMI)", "Number of vehicles in household", "Cheese intake", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]
    for chr in 1:3
        dataset = DataFrame(Arrow.Table(joinpath("results", "datasets", "subset_chr$(chr).data.arrow")))
        @test nrow(dataset) == 500
        columnnames = names(dataset)
        @test issubset(traits_and_pcs, columnnames)
        variant_columns = setdiff(columnnames, traits_and_pcs)
        # All variants should be from the specified chromosome
        @test all(startswith(colname, string(chr)) for colname in variant_columns)
    end

    # Check QQ plot exists
    @test isfile(joinpath("results", "QQ.png"))

    # Check pve output files
    @test isfile(joinpath("results", "pve", "pve.subset_chr1.txt"))
    @test isfile(joinpath("results", "pve", "pve.subset_chr2.txt"))
    @test isfile(joinpath("results", "pve", "pve.subset_chr3.txt"))

    # Check results - should only contain the 9 variants from subset_bed_files.txt
    results = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    @test size(results, 1) > 0
    
    # Collect all unique variants from results
    variants = Set{Symbol}()
    failed_results = []
    estimator_names = filter(x -> !endswith(x, "PVALUE"), names(results))
    for estimator_name in estimator_names
        for Ψ̂ in results[!, estimator_name]
            if Ψ̂ isa TMLECLI.FailedEstimate
                push!(failed_results, Ψ̂)
            end
            Ψ = Ψ̂.estimand
            @test Ψ isa JointEstimand
            @test get_outcome(Ψ) == Symbol("Body mass index (BMI)")
            variant = only(TargeneCore.get_treatments(Ψ))
            push!(variants, variant)
        end
    end
    
    # Check that we only have variants from the subset file (9 variants across 3 chromosomes)
    @test length(variants) <= 9
    @test Set(string(v)[1] for v ∈ variants) == Set(['1', '2', '3'])
    @test size(failed_results, 1) == 0

    # Check summary file
    summary_results = YAML.load_file(joinpath("results", "results.summary.yaml"))
    @test size(summary_results, 1) == size(results, 1)

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 100
end

end

true
