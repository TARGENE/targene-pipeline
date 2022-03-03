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
    
    grm_ids = UKBBEpistasisPipeline.GRMIDs(joinpath("data", "grm", "test.grm.id"))

    build_result_file(grm_ids)

    results_file = jldopen("results.hdf5", "a+")

    phenotypes, inf_curves, n_observations = UKBBEpistasisPipeline.build_influence_curves(results_file, grm_ids)
    @test phenotypes == ["cancer", "bmi"]
    @test size(inf_curves) == (size(grm_ids, 1), 2, 2)
    @test n_observations == [184, 174]


    # The operation is converting to Float32
    # All missing elements are at the enf of the grm ids
    qrs = results_file["cancer"]["queryreports"]
    @test expected_inf_curve(qrs, 1, 10) == inf_curves[:, 1, 1]
    @test expected_inf_curve(qrs, 2, 10) == inf_curves[:, 2, 1]

    qrs = results_file["bmi"]["queryreports"]
    @test expected_inf_curve(qrs, 1, 20) == inf_curves[:, 1, 2]
    @test expected_inf_curve(qrs, 2, 20) == inf_curves[:, 2, 2]

    close(results_file)
    rm("results.hdf5")
end


@testset "Test bit_distance" begin
    grm_sample = Float32[-0.6, -0.4, -0.25, -0.3, -0.1, 0.1, 0.3, 0.5, 0.2, 0.6]
    nτs = 5
    τs = UKBBEpistasisPipeline.default_τs(nτs)
    @test τs ≈ Float32[2.0  1.0  0.666667  0.5  0.4]
    d = UKBBEpistasisPipeline.bit_distances(grm_sample, τs)
    @test d == [1  0  0  0  0
                1  0  0  0  0
                1  0  0  0  0
                1  0  0  0  0
                1  0  0  0  0
                1  1  0  0  0
                1  1  1  1  1
                1  1  1  1  1
                1  1  1  0  0
                1  1  1  1  1]
end


@testset "Test compute_variances" begin
    n_phenotypes = 3
    n_queries = 2
    n_samples = 5
    nτs = 4
    n_obs = [3, 5, 4]
    τs = UKBBEpistasisPipeline.default_τs(nτs)
    # The GRM has 15 lower triangular elements
    grm = [0.4, 0.1, 0.5, 0.2, -0.2, 0.6, 0.3, -0.6, 
            0.4, 0.3, 0.6, 0.3, 0.7, 0.3, 0.1]
    influence_curves = zeros(Float32, n_samples, n_queries, n_phenotypes)
    influence_curves[:, 1, 1] = [0.1, 0., 0.1, 0.3, 0.]
    influence_curves[:, 1, 2] = [0.1, 0.2, 0.1, 0.3, 0.2]
    influence_curves[:, 1, 3] = [0.1, 0., 0.1, 0.3, 0.2]
    influence_curves[:, 2, 1] = [0.1, 0., 0.1, 0.3, 0.]
    influence_curves[:, 2, 2] = [0.1, 0.2, 0.1, 0.3, 0.2]
    influence_curves[:, 2, 3] = [0.1, 0., 0.1, 0.3, 0.2]
                  
    
    @enter variances = UKBBEpistasisPipeline.compute_variances(influence_curves, grm, τs, n_obs)
    @test size(variances) == (nτs, n_queries, n_phenotypes)

    # when τ=2, all elements are used
    for phenotype_idx in 1:n_phenotypes
        for query_idx in n_queries
            s = sum(influence_curves[:, query_idx, phenotype_idx])
            var = sum(s*influence_curves[i, query_idx, phenotype_idx] for i in 1:n_samples)/n_obs[phenotype_idx]
            @test variances[1, query_idx, phenotype_idx] ≈ var
        end
    end

    # Decreasing variances with τ as all inf curves are positives
    for nτ in 2:nτs
        @test all(variances[nτ, :, :] .<= variances[nτ-1, :, :])
    end

end


@testset "Test sieve_variance_plateau" begin

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