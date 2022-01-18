module TestGRM

using Test
using UKBBEpistasisPipeline
using SnpArrays
using Serialization

@testset "Test computeGRM" begin
    # As I use multiple times the same file, the GRM constructed from those files should 
    # be the same as the GRM computed from only 1 file
    outdir = "grmmatrix"
    parsed_args = Dict(
        "outdir" => outdir,
        "bed-files" => [SnpArrays.datadir("mouse.bed"), SnpArrays.datadir("mouse.bed")]
    )

    snparray = SnpArray(SnpArrays.datadir("mouse.bed"))

    # Taken from the source of the package but this uses matrix multiplication
    basegrm = grm(snparray)

    computeGRM(parsed_args)
    # Check proper amount of lines are written
    @test length(readdir(outdir)) == 1940
    # Check a few lines at random
    for i in (1, 20, 365, 1856)
        grm_line = deserialize(joinpath(outdir, "grm_$i.jls"))
        @test grm_line ≈ basegrm[i, :]
    end
    rm(outdir, recursive=true)
end

end

true