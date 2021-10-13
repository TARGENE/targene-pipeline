module TestGenerateQueries

using Test
using UKBBEpistasisPipeline
using DataFrames
using TOML

chrpath(x) = joinpath("data", "ukbb", "imputed" ,"ukbb_chr$x.bgen")


@testset "Test get_minor_major" begin
    parsed_args = Dict(
        "sample-chrom-file" => chrpath(8),
        "trans-actors-file" => joinpath("data", "trans_actors_fake.csv"),
        "asb-file" => joinpath("data", "filtered_asb_fake.csv"),
        "exclude" => joinpath("data", "pb_snps.csv")
        )
    trans_actors, asbs = UKBBEpistasisPipeline.get_minor_major(parsed_args)

    @test trans_actors == DataFrame(
        ID        = ["RSID_102", "RSID_2"],
        CHROM     = ["chr12", "chr12"],
        MAJOR     = ["A", "G"],
        MINOR     = ["G", "A"],
        CHROMPATH = ["data/ukbb/imputed/ukbb_chr12.bgen", "data/ukbb/imputed/ukbb_chr12.bgen"]
    )

    @test asbs == DataFrame(
        ID        = ["RSID_198", "RSID_99"],
        CHROM     = ["chr12", "chr12"],
        MAJOR     = ["G", "G"],
        MINOR     = ["A", "A"],
        CHROMPATH = ["data/ukbb/imputed/ukbb_chr12.bgen", "data/ukbb/imputed/ukbb_chr12.bgen"]
    )
end


@testset "Test generate_queries" begin
    outdir = joinpath("data", "queries")

    mkdir(outdir)

    parsed_args = Dict(
        "mode" => "frequency",
        "out" => outdir,
        "threshold" => 0.8,
        "sample-chrom-file" => chrpath(8),
        "trans-actors-file" => joinpath("data", "trans_actors_fake.csv"),
        "asb-file" => joinpath("data", "filtered_asb_fake.csv"),
        "exclude" => joinpath("data", "pb_snps.csv")
        )

    generate_queries(parsed_args)

    # Check generated filenames
    query_files = readdir(outdir)
    @test query_files == ["query_RSID_198_RSID_102.toml",
                          "query_RSID_198_RSID_2.toml",
                          "query_RSID_99_RSID_102.toml",
                          "query_RSID_99_RSID_2.toml"]

    # Check content of the files
    expected_query = Dict(
        "RSID_2" => "GG -> GA",
        "RSID_198" => "GG -> GA",
        "RSID_102" => "AA -> AG",
        "RSID_99" => "GG -> GA"
    )
    for file in query_files
        queries = TOML.parsefile(joinpath(outdir, file))

        @test queries["threshold"] == 0.8

        snps = queries["SNPS"]
        @test length(snps) == 2

        for (snp, path) in snps
            @test path == chrpath(12)
        end

        query = queries["QUERY"]
        @test length(query) == 2
        for (snp, q) in query
            @test q == expected_query[snp]
        end
    end

    # Cleanup
    rm(outdir; recursive=true)
end

end;

true