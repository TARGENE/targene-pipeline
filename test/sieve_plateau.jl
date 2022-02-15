module TestSievePlateau

using Test
using UKBBEpistasisPipeline
using DataFrames
using CSV 
using JLD2
using TMLE


@testset "Test build_influence_curves" begin
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

@testset "Test generate_indices" begin
    n_samples = 10
    chunk_sizes = [14, 22, 19]
    gen = UKBBEpistasisPipeline.generate_indices(chunk_sizes, n_samples)
    @test take!(gen) == (5, 4)
    @test take!(gen) == (8, 8)
    @test take!(gen) == (10, 10)
    @test_throws TaskFailedException take!(gen)
end

@testset "Test chunk_pairwise_inf_curve" begin
    inf_curve = collect(1:10)
    first_i, first_j = 1, 1
    last_i, last_j = 3, 1
    pairwise_ic = UKBBEpistasisPipeline.chunk_pairwise_inf_curve(inf_curve, first_i, first_j, last_i, last_j)
    @test pairwise_ic == [1., 2., 4., 3.]
    
    first_i, first_j = 3, 2
    last_i, last_j = 6, 6
    pairwise_ic = UKBBEpistasisPipeline.chunk_pairwise_inf_curve(inf_curve, first_i, first_j, last_i, last_j)
    @test pairwise_ic == 
        [6., 9., 4., 8., 12., 16., 5., 10., 15., 20., 25., 6.]

end

@testset "Test" begin
    grm_chunk_file = "GRM_1.bin"
    open(grm_chunk_file, "w") do io 
        write(io, 10)
        write(io, Float32[-1.1, -0.8, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.8, 1.1])
    end
    nτs = 5
    # τs = reshape([2/i for i in 1:nτs], 1, nτs)
    d = UKBBEpistasisPipeline.bit_distance(grm_chunk_file, nτs)
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
    rm(grm_chunk_file)
end


@testset "Test compute_variances" begin
    # Sizes
    n_samples = 10
    n_queries = 2
    n_phenotypes = 3
    # Generate GRM and split
    grm_size = n_samples*(n_samples+1) ÷ 2
    h = 2/(grm_size-1)
    grm = collect(-1:h:1)

    previous = 1
    for (id, number) in [(1, 14), (2, 22), (3, 19)]
        open("GRM_$id.bin", "w") do io 
            write(io, number)
            write(io, grm[previous:previous+number-1])
        end
        previous += number
    end

    # Compute variances
    grmfiles = ["GRM_1.bin", "GRM_2.bin", "GRM_3.bin"]
    influence_curves = reshape(collect(1:n_phenotypes*n_samples*n_queries), 
                        n_samples, 
                        n_phenotypes, 
                        n_queries)
    nτs = 5
    variances = UKBBEpistasisPipeline.compute_variances(influence_curves, nτs, grmfiles)
    @test size(variances) == (n_phenotypes, n_queries, nτs)

    # For the first τ (max distance==2), all elements are taken into account
    for phenotype_idx in 1:n_phenotypes
        for query_idx in 1:n_queries
            inf_phen_1 = Float32.(influence_curves[:, phenotype_idx, query_idx])    
            @test variances[phenotype_idx, query_idx, 1] == 
                sum(inf_phen_1[i]*inf_phen_1[j] 
                        for i in 1:n_samples for j in 1:i)

        end
    end

    # Decreasing variances with τ as all inf curves are positives
    for nτ in 2:nτs
        @test all(variances[:, :, nτ] .<= variances[:, :, nτ-1])
    end

    # Clean
    for id in (1,2,3)
        rm("GRM_$id.bin")
    end
end

@testset "Test compute_variances" begin
    parsed_args = Dict(
        "grm-ids" => "data/grm/grm_bis.id",
        "results" => "data/RSID_10_RSID_100.hdf5",
        "nτs" => 10
    )
end

end