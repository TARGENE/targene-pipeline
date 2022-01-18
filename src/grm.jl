
using Serialization 

function update(outdir, i, grmᵢ, nbSNPs, finalize)
    filepath = joinpath(outdir, "grm_$i.jls")
    if isfile(filepath)
        old_grmᵢ = deserialize(filepath)
        grmᵢ = grmᵢ + old_grmᵢ
    end
    grmᵢ = finalize ? grmᵢ / (2*(nbSNPs-1)) : grmᵢ
    serialize(filepath, grmᵢ)
end


function computeGRM(parsed_args)
    outdir = parsed_args["outdir"]
    bedfiles = parsed_args["bed-files"]
    nfiles = length(bedfiles)
    nbSNPs = 0
    mkdir(outdir)
    for (fileid, bedfile) in enumerate(bedfiles)
        println("Aggregating GRM from file: $bedfile")
        snparray = SnpArray(bedfile)
        n, p = size(snparray)
        # Update total amount of SNPs for normalization
        nbSNPs += p

        # Compute centered/normalized genotype counts
        G = Matrix{Float32}(undef, n, p)
        Base.copyto!(G, snparray, model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
        
        # Update individual GRMᵢ in multithreading mode
        Threads.@threads for i in 1:n
            grmᵢ = G*G[i, :]
            update(outdir, i, grmᵢ, nbSNPs, fileid==nfiles)
        end
    end
end