module TestConfoundersPreparation

using Test
using SnpArrays
using UKBBEpistasisPipeline
using DataFrames
using CSV

@testset "Various functions" begin
    # Test issnp
    @test UKBBEpistasisPipeline.issnp("A") == true
    @test UKBBEpistasisPipeline.issnp("AA") == false
    
    #Test bounds
    @test UKBBEpistasisPipeline.bounds(0, -10, 10) == (-10^4, 10^4)
    @test UKBBEpistasisPipeline.bounds(15000, 1000, 30000) == (1000, 30000)

    # Test notin_ldblocks
    ldblocks = DataFrame(chr        = [1, 3],
                        lower_bound = [10, 19],
                        upper_bound = [100, 109])
    ldblocks = groupby(ldblocks, :chr)

    snp_info = DataFrame(chromosome=[1, 1, 2, 3, 3],
                        position=[21, 5, 9, 102, 4])
                        \
    expected_snps = DataFrame(chromosome = [1, 2, 3],
                              position   = [5, 9, 4])
    @test filter(x -> UKBBEpistasisPipeline.notin_ldblocks(x, ldblocks), snp_info) == expected_snps
end

@testset "Test filter_chromosome" begin
    parsed_args = Dict(
        "input"  => SnpArrays.datadir("mouse"),
        "output" => joinpath("data", "filtered-mouse"),
        "qcfile" => joinpath("data", "ukbb", "qcfile.txt"),
        "ld-blocks" => joinpath("data", "VDR_LD_blocks.txt"),
        "maf-threshold" => 0.31
    )
    filter_chromosome(parsed_args)

    filtered = SnpData(parsed_args["output"])

    @test filtered.snp_info.snpid == ["rs13476318"]
    @test size(filtered.snparray) == (1940, 1)

    # Clean
    for ext in [".bed", ".bim", ".fam"]
        rm(parsed_args["output"]*ext)
    end
end

@testset "Test merge_beds" begin
    genotypes_dir = joinpath("data", "ukbb", "genotypes")
    chr_prefix = joinpath(genotypes_dir, "ukbb")

    parsed_args = Dict(
            "input"  => chr_prefix,
            "output" => joinpath(genotypes_dir, "ukbb_merged")
        )
    merge_beds(parsed_args)
    merged = SnpData(parsed_args["output"])

    @test length(unique(merged.snp_info.chromosome)) == 3
    # Clean

    for ext in [".bed", ".bim", ".fam"]
        rm(parsed_args["output"]*ext)
    end

end


@testset "Test adapt_flashpca" begin
    parsed_args = Dict(
        "input" => joinpath("data", "flashpca_output.txt"),
        "output" => "flashpca_after_adapt.csv",
    )
    adapt_flashpca(parsed_args)

    output = CSV.File(parsed_args["output"]) |> DataFrame

    @test names(output) == ["SAMPLE_ID",
                            "PC1",
                            "PC2",
                            "PC3",
                            "PC4",
                            "PC5",
                            "PC6"]
    @test size(output) == (7, 7)
    
    rm(parsed_args["output"])
end


end;

true