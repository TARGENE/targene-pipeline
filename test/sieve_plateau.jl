module TestSievePlateau

using Test
using UKBBEpistasisPipeline
using DataFrames
using CSV 
using JLD2
using TMLE

function build_result_file(grm_ids;path="results.hdf5")
    # The queries are not important here
    query = TMLE.Query(case=(t₁="CG", t₂="TT"), control=(t₁="GG", t₂="TA"), name="MyQuery")
    jldopen(path, "w") do file
        # First phenotype
        cancer = JLD2.Group(file, "cancer")
        n_samples = size(grm_ids, 1) - 10
        cancer["sample_ids"] = grm_ids.SAMPLE_ID[1:n_samples]
        cancer["queryreports"] = [TMLE.QueryReport(query, collect(1:n_samples), 0.1, 0.2),
                                  TMLE.QueryReport(query, collect(n_samples:2n_samples-1), 0.1, 0.2)]
        # Second phenotype
        bmi = JLD2.Group(file, "bmi")
        n_samples = size(grm_ids, 1) - 20
        bmi["sample_ids"] = grm_ids.SAMPLE_ID[1:n_samples]
        bmi["queryreports"] = [TMLE.QueryReport(query, collect(1:n_samples), 0.1, 0.2),
                                TMLE.QueryReport(query, collect(n_samples:2n_samples-1), 0.1, 0.2)]

    end
end


@testset "Test build_influence_curves" begin
    expected_inf_curve(queryreports, index, nbzeros) = 
        convert(Vector{Float32}, vcat(queryreports[index].influence_curve, zeros(nbzeros)))
    
    grm_ids = UKBBEpistasisPipeline.load_grm_ids("data/grm/GRM.ids.csv")

    build_result_file(grm_ids)

    results_file = jldopen("results.hdf5", "a+")

    phenotypes, inf_curves = UKBBEpistasisPipeline.build_influence_curves(results_file, grm_ids)
    @test phenotypes == ["cancer", "bmi"]
    @test size(inf_curves) == (size(grm_ids, 1), 2, 2)

    # The operation is converting to Float32
    # All missing elements are at the enf of the grm ids
    qrs = results_file["cancer"]["queryreports"]
    @test expected_inf_curve(qrs, 1, 10) == inf_curves[:, 1, 1]
    @test expected_inf_curve(qrs, 2, 10) == inf_curves[:, 1, 2]

    qrs = results_file["bmi"]["queryreports"]
    @test expected_inf_curve(qrs, 1, 20) == inf_curves[:, 2, 1]
    @test expected_inf_curve(qrs, 2, 20) == inf_curves[:, 2, 2]

    close(results_file)
    rm("results.hdf5")
end


@testset "Test bit_distance" begin
    grm_sample = Float32[-1.1, -0.8, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.8, 1.1]
    nτs = 5
    # τs = reshape([2/i for i in 1:nτs], 1, nτs)
    d = UKBBEpistasisPipeline.bit_distances(grm_sample, nτs)
    @test d == [1  0  0  0  0
                1  0  0  0  0
                1  0  0  0  0
                1  0  0  0  0
                1  0  0  0  0
                1  1  0  0  0
                1  1  0  0  0
                1  1  1  1  0
                1  1  1  1  1
                1  1  1  1  1]
end


@testset "Test compute_variances" begin
    # Compute variances
    grm, grm_ids = UKBBEpistasisPipeline.readGRM("data/grm/GRM")
    n_phenotypes = 2
    n_queries = 2
    n_samples = size(grm_ids, 1)
    influence_curves = reshape(collect(1:n_phenotypes*n_samples*n_queries), 
                        n_samples, 
                        n_phenotypes, 
                        n_queries)
    influence_curves = convert(Array{Float32, 3}, influence_curves)                 
    nτs = 5
    variances = UKBBEpistasisPipeline.compute_variances(influence_curves, grm, nτs)
    @test size(variances) == (n_phenotypes, n_queries, nτs)

    # For the first τ (max distance==2), all elements are taken into account
    for phenotype_idx in 1:n_phenotypes
        for query_idx in 1:n_queries
            inf_phen_1 = influence_curves[:, phenotype_idx, query_idx]  
            @test variances[phenotype_idx, query_idx, 1] ≈ 
                sum(inf_phen_1[i]*inf_phen_1[j] 
                        for i in 1:n_samples for j in 1:n_samples)

        end
    end

    # Decreasing variances with τ as all inf curves are positives
    for nτ in 2:nτs
        @test all(variances[:, :, nτ] .<= variances[:, :, nτ-1])
    end

end

@testset "Test compute_variances" begin

    grm_ids = UKBBEpistasisPipeline.load_grm_ids("data/grm/GRM.ids.csv")
    build_result_file(grm_ids;path="results.hdf5")
    parsed_args = Dict(
        "grm-prefix" => "data/grm/GRM",
        "results" => "results.hdf5",
        "nb-estimators" => 10
    )

    sieve_variance_plateau(parsed_args)

    results_file = jldopen("results.hdf5", "a+")

    @test size(results_file["cancer"]["variances"]) == (2, 10)
    @test size(results_file["bmi"]["variances"]) == (2, 10)

    rm("results.hdf5")
end

end