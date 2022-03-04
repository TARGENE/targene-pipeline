module TestSievePlateau

using Test
using UKBBEpistasisPipeline
using DataFrames
using CSV 
using JLD2
using TMLE
using CategoricalArrays
using MLJLinearModels
using MLJBase
using Random

function build_result_file(grm_ids; path="results_test.hdf5")
    rng = Xoshiro(0)
    n = size(grm_ids, 1)
    T = (t₁=categorical(rand(rng, [0, 1], n)),)
    W = (w₁=rand(rng, n), w₂=rand(rng, n))
    y₁ = rand(rng, n)
    y₂ = 2convert(Vector{Float64}, T.t₁) + 0.2rand(rng, n)
    query_1 = Query(case=(t₁=1,), control=(t₁=0,))
    query_2 = Query(case=(t₁=0,), control=(t₁=1,))
    tmle = TMLEstimator(
        LinearRegressor(), 
        LogisticClassifier(),
        query_1,
        query_2)
    mach₁ = machine(tmle, T, W, y₁)
    fit!(mach₁, verbosity=0)
    mach₂ = machine(tmle, T, W, y₂)
    fit!(mach₂, verbosity=0)

    jldopen(path, "w") do io
        cancer = JLD2.Group(io, "cancer")
        cancer["queryreports"] = getqueryreports(mach₁)
        cancer["sample_ids"] = string.(grm_ids.SAMPLE_ID)

        bmi = JLD2.Group(io, "bmi")
        bmi["queryreports"] = getqueryreports(mach₂)
        bmi["sample_ids"] = string.(grm_ids.SAMPLE_ID)
    end
end

@testset "Test build_work_list" begin
    grm_ids = UKBBEpistasisPipeline.GRMIDs("data/grm/test.grm.id")
    path = "results_test.hdf5"
    build_result_file(grm_ids; path=path)
    results_file = jldopen(path, "a")
    influence_curves, n_obs, phenotype_query_pairs = UKBBEpistasisPipeline.build_work_list(results_file, grm_ids; pval=0.05)
    bmi_qrs = results_file["bmi"]["queryreports"]
    @test influence_curves == convert(Matrix{Float32}, [bmi_qrs[1].influence_curve'; bmi_qrs[2].influence_curve'])
    @test n_obs == [194, 194]
    @test phenotype_query_pairs == ["bmi"=> 1, "bmi"=>2]
    rm(path)
end

@testset "Test bit_distance" begin
    sample_grm = Float32[-0.6, -0.4, -0.25, -0.3, -0.1, 0.1, 0.3, 0.5, 0.2, 0.6]
    nτs = 5
    τs = UKBBEpistasisPipeline.default_τs(nτs)
    @test τs ≈ Float32[2.0, 1.0,  0.666667,  0.5,  0.4]
    d = UKBBEpistasisPipeline.bit_distances(sample_grm, τs)
    @test d == [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
                0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  1.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  1.0
                0.0  0.0  0.0  0.0  0.0  0.0  1.0  1.0  0.0  1.0]
end

@testset "Test aggregate_variances" begin
    # 2 influence curves containing 5 individuals
    influence_curves = [1. 2. 3. 4. 5.
                        6. 7. 8. 9. 10.]
    # distance indicator with 3 τs and corresponding to row 4
    indicator = [1. 0. 0. 1.
                 0. 0. 1. 1.
                 1. 0. 1. 1.]
    sample = 4
    var = UKBBEpistasisPipeline.aggregate_variances(influence_curves, indicator, sample)
    @test var == [24.0  189.0
                  40.0  225.0
                  48.0  333.0]
end

@testset "Test normalize!" begin
    # 2 τs and 3 curves
    n_obs = [10, 10, 100]
    variances = [1. 2. 3.
                 4. 5. 6.]
    UKBBEpistasisPipeline.normalize!(variances, n_obs)
    @test variances == [0.1 0.2 0.03
                        0.4 0.5 0.06]
end

@testset "Test compute_variances" begin
    n_curves = 3
    n_samples = 5
    nτs = 4
    n_obs = [3, 4, 4]
    τs = UKBBEpistasisPipeline.default_τs(nτs)
    # The GRM has 15 lower triangular elements
    grm = [0.4, 0.1, 0.5, 0.2, -0.2, 0.6, 0.3, -0.6, 
            0.4, 0.3, 0.6, 0.3, 0.7, 0.3, 0.1]
    influence_curves = [0.1 0. 0.1 0.3 0.
                        0.1 0.2 0.1 0.3 0.2
                        0.1 0. 0.1 0.3 0.2]
                  
    
    variances = UKBBEpistasisPipeline.compute_variances(influence_curves, grm, τs, n_obs)
    @test size(variances) == (nτs, n_curves)

    # when τ=2, all elements are used
    for curve_id in 1:n_curves
        s = sum(influence_curves[curve_id, :])
        var = sum(s*influence_curves[curve_id, i] for i in 1:n_samples)/n_obs[curve_id]
        @test variances[1, curve_id] ≈ var
    end

    # Decreasing variances with τ as all inf curves are positives
    for nτ in 2:nτs
        @test all(variances[nτ, :] .<= variances[nτ-1, :])
    end

end

@testset "Test sieve_variance_plateau" begin
    grm_ids = UKBBEpistasisPipeline.GRMIDs("data/grm/test.grm.id")
    path = "results_test.hdf5"
    build_result_file(grm_ids; path=path)
    parsed_args = Dict(
        "grm-prefix" => "data/grm/test.grm",
        "results" => path,
        "nb-estimators" => 10
    )

    sieve_variance_plateau(parsed_args)

    results_file = jldopen(path, "a+")

    @test haskey(results_file["cancer"], "sieve_variances") == false
    @test size(results_file["bmi"]["sieve_variances"]["1"]) == (10,)
    @test size(results_file["bmi"]["sieve_variances"]["2"]) == (10,)

    rm(path)
end


@testset "Test grm_rows_bounds" begin
    n_samples = 5
    grm_bounds = UKBBEpistasisPipeline.grm_rows_bounds(n_samples)
    @test grm_bounds == [1 => 1
                         2 => 3
                         4 => 6
                         7 => 10
                         11 => 15]
end

# @testset "Test perf" begin
#     sample = 200_000
#     n_curves = 100
#     nτs = 100
#     τs = UKBBEpistasisPipeline.default_τs(nτs)
#     sample_grm = rand(Float32, sample)
#     indicator = UKBBEpistasisPipeline.bit_distances(sample_grm, τs)
#     D = rand(Float32, n_curves, 2sample)
#     @btime UKBBEpistasisPipeline.aggregate_variances(D, indicator, sample)

# end

end