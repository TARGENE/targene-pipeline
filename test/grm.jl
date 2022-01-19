module TestGRM

using Test
using UKBBEpistasisPipeline
using SnpArrays
using HDF5


@testset "Test computeGRM" begin
    # As I use multiple times the same file, the GRM constructed from those files should 
    # be the same as the GRM computed from only 1 file
    outfile = "grmmatrix.hdf5"
    parsed_args = Dict(
        "outfile" => outfile,
        "bed-files" => [SnpArrays.datadir("mouse.bed"), SnpArrays.datadir("mouse.bed")]
    )

    snparray = SnpArray(SnpArrays.datadir("mouse.bed"))

    # Taken from the source of the package but this uses matrix multiplication
    basegrm = grm(snparray)

    computeGRM(parsed_args)
    # Check proper amount of lines are written
    file = h5open(outfile)

    @test file["GRM"][:,:] ≈ basegrm

    rm(outfile)
end

end

true
