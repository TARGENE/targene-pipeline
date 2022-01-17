
using Serialization 

function update(outdir, i, grmᵢ, R, finalize)
    filepath = joinpath(outdir, "grm_$i.jls")
    if isfile(filepath)
        old_grmᵢ = deserialize(filepath)
        grmᵢ = grmᵢ + old_grmᵢ
    end
    grmᵢ = finalize ? grmᵢ / (R-1) : grmᵢ
    serialize(filepath, grmᵢ)
end


function computeGRM(outdir, bedfiles)
    nfiles = length(bedfiles)
    R = 0
    mkdir(outdir)
    for (fileid, bedfile) in enumerate(bedfiles)
        snparray = SnpArray(bedfile)
        n, p = size(snparray)

        # Update total amount of SNPs for normalization
        R += p

        # Compute centered genotype counts
        allelecounts = counts(snparray, dims=1)
        a2freqs = (2allelecounts[end, :] + allelecounts[end-1, :]) / 2n
        a2copies = convert(Matrix{Float32}, snparray) .- 2repeat(transpose(a2freqs), n)
        replace!(a2copies, NaN => 0)
        den = 2a2freqs.*(1 .- a2freqs)
        
        # Update individual GRMᵢ in multithreading mode
        for i in 1:n
            grmᵢ = a2copies*(a2copies[i, :] ./ den)
            update(outdir, i, grmᵢ, R, fileid==nfiles)
        end
    end
end