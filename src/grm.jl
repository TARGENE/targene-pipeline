
function computeGRM(parsed_args)
    bedfiles = parsed_args["bed-files"]
    nfiles = length(bedfiles)

    nbSNPs = 0
    file = h5open(parsed_args["outfile"], "w")
    GRM = nothing

    for (fileid, bedfile) in enumerate(bedfiles)
        println("Aggregating GRM from file: $bedfile")
        snparray = SnpArray(bedfile)
        n, p = size(snparray)
        if GRM === nothing
            GRM = create_dataset(file, "GRM", datatype(Float32), dataspace(n, n), chunk=(n, 1), blosc=3)
        end
        # Update total amount of SNPs for normalization
        nbSNPs += p

        # Compute centered/normalized genotype counts
        G = Matrix{Float32}(undef, n, p)
        Base.copyto!(G, snparray, model=ADDITIVE_MODEL, impute=true, center=true, scale=true)
        
        # Update individual GRMᵢ
        for i in 1:n
            GRM[:, i] += G*G[i, :]
            # Normalize after last file
            if (fileid == nfiles)
                GRM[:, i] /= (2*(nbSNPs-1))
            end
        end
    end
    close(file)
end


