module TestGWASReportBinary

using Test
using DataFrames
using CSV
using TMLE
using Arrow
using JLD2
using TMLECLI
using TargeneCore
using YAML

args = length(ARGS) > 0 ? ARGS : ["-profile", "local", "-resume"]

# Exercises the TarGWAS report on a 2-SNP GWAS of a binary UKB trait
# (ICD-10 I10, essential hypertension). Beyond the usual estimation checks,
# this asserts the report artefacts: the HTML page, the enriched per-estimator
# CSVs (carrying genotype counts + BIM-corrected chrom/pos) and the LD matrix.
@testset "Test ukb_gwas_report_binary.config" begin
    cmd = `nextflow run main.nf -c test/configs/ukb_gwas_report_binary.config $args`
    @info string("The following command will be run:\n", cmd)

    r = run(cmd)
    @test r.exitcode == 0

    # Outcome is the binary I10 trait; the two extra traits are covariates.
    outcome = Symbol("I10 Essential hypertension")
    traits_and_pcs = ["SAMPLE_ID", "I10 Essential hypertension",
        "Number of vehicles in household", "Cheese intake",
        "PC1", "PC2", "PC3", "PC4", "PC5", "PC6"]

    # Only chr1 is subset to two SNPs for the analysis dataset.
    dataset = DataFrame(Arrow.Table(joinpath("results", "datasets", "subset.ukb_chr1.data.arrow")))
    @test nrow(dataset) == 500
    columnnames = names(dataset)
    @test issubset(traits_and_pcs, columnnames)
    variant_columns = setdiff(columnnames, traits_and_pcs)
    @test length(variant_columns) == 2
    @test all(startswith(colname, "1") for colname in variant_columns)

    # Estimation results: only the two requested SNPs.
    results = jldopen(io -> io["results"], joinpath("results", "results.hdf5"))
    @test size(results, 1) > 0
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
            @test get_outcome(Ψ) == outcome
            push!(variants, only(TargeneCore.get_treatments(Ψ)))
        end
    end
    @test length(variants) == 2
    @test all(startswith(string(v), "1:") for v in variants)
    @test size(failed_results, 1) == 0

    # ----------------------------------------------------------------------
    # TarGWAS report artefacts.
    # ----------------------------------------------------------------------
    report_dir = joinpath("results", "report")
    @test isfile(joinpath(report_dir, "report.html"))

    html = read(joinpath(report_dir, "report.html"), String)
    @test occursin("TarGene GWAS report", html)
    @test occursin("plotly", lowercase(html))
    # The interactive hits table + genotype-count threshold filter must render.
    @test occursin("hits-table", html)
    @test occursin("tbl-threshold", html)
    @test occursin("\"has_geno\":true", html)

    # Enriched per-estimator CSVs carry genotype-count columns.
    tabulated_dirs = filter(p -> isdir(p) && endswith(p, "_tabulated"),
        readdir(report_dir; join=true))
    @test !isempty(tabulated_dirs)
    tabulated = first(tabulated_dirs)
    csvs = filter(f -> startswith(basename(f), "GWAS_results_") && endswith(f, ".csv"),
        readdir(tabulated; join=true))
    @test !isempty(csvs)
    csv = DataFrame(CSV.File(first(csvs); delim='\t'))
    @test nrow(csv) == 2                       # two SNPs
    @test "case_geno_counts" in names(csv)
    @test "ctrl_geno_counts" in names(csv)
    @test all(occursin(":", s) for s in skipmissing(csv.case_geno_counts))
    # chrom/pos resolved from the BIM (chr1), not zeroed.
    @test all(csv.chrom .== 1)
    @test all(csv.pos .> 0)

    # LD matrix is persisted for downstream analysis.
    @test isfile(joinpath(tabulated, "LD_matrix.csv"))

    # Check properly resumed
    resume_time = @elapsed run(cmd)
    @test resume_time < 100
end

end

true
