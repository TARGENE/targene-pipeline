module TestSievePlateau

using Test
using UKBBEpistasisPipeline
using DataFrames
using CSV 
using JLD2
using TMLE


@testset "Test influence curves re-build" begin
    grm_ids = UKBBEpistasisPipeline.load_grm_ids("data/grm/grm_bis.id")
    results_file = jldopen("data/RSID_10_RSID_100.hdf5", "a+")

    phenotypes, inf_curves = UKBBEpistasisPipeline.build_influence_curves(results_file, grm_ids)
    @test phenotypes == ["categorical_phenotype", "continuous_phenotype"]
    @test size(inf_curves) == (size(grm_ids, 1), 2, 2)

    # The operation is converting to Float32
    # All missing elements are at the enf of the grm ids
    qrs = results_file["categorical_phenotype"]["queryreports"]
    @test convert(Vector{Float32}, qrs[1].influence_curve) == collect(skipmissing(inf_curves[:, 1, 1]))
    @test convert(Vector{Float32}, qrs[2].influence_curve) == collect(skipmissing(inf_curves[:, 1, 2]))

    qrs = results_file["continuous_phenotype"]["queryreports"]
    @test convert(Vector{Float32}, qrs[1].influence_curve) == collect(skipmissing(inf_curves[:, 2, 1]))
    @test convert(Vector{Float32}, qrs[2].influence_curve) == collect(skipmissing(inf_curves[:, 2, 2]))

    close(results_file)
end

@testset "Test variances computations" begin
    parsed_args = Dict(
        "grm-ids" => "data/grm/grm_bis.id",
        "results" => "data/RSID_10_RSID_100.hdf5",
        "nτs" => 10
    )
end

end