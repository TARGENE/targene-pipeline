module TestGRM

using Test
using UKBBEpistasisPipeline
using SnpArrays

@testset "Test computeGRM" begin
    outdir = "grmmatrix"
    bedfiles = [SnpArrays.datadir("mouse.bed")]
    snparray = SnpArray(SnpArrays.datadir("mouse.bed"))
    # Taken from the source of the package
    m, n = size(snparray)
    G = Matrix{Float32}(undef, m, n)
    Base.copyto!(G, snparray, model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
    Φ = G * transpose(G)
    Φ ./= 2n

    basegrm = grm(snparray)

    computeGRM(outdir, bedfiles)

    grm₁ = deserialize(joinpath(outdir, "grm_1.jls"))
    grm₁/2
    rm(outdir, recursive=true)
end

end

true