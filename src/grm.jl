
function grm_from_bedfile(bedfile)
    println("Computing partial GRM from file: $bedfile")
    snparray = SnpArray(bedfile)
    n, p = size(snparray)
    grm(snparray), (n, p)
end

updateGRM!(grm, partial_grm, p) = grm[:, :] += partial_grm*p

function initialize_grm(grmfile, bedfile)
    initial_grm, (n, p) = grm_from_bedfile(bedfile)
    
    grm = create_dataset(grmfile, "GRM", datatype(Float32), dataspace(n, n), chunk=(n, 1), blosc=3)
    updateGRM!(grm, initial_grm, p)

    return grm, p
end


function computeGRM(parsed_args)
    outfile = h5open(parsed_args["outfile"], "w")
    bedfiles = parsed_args["bed-files"]

    # Initialize the GRM from the first bedfile
    grm, nbSNPs = initialize_grm(outfile, bedfiles[1])

    for bedfile in bedfiles[2:end]
        partial_grm, (n, p) = grm_from_bedfile(bedfile)
        updateGRM!(grm, partial_grm, p)
        # Update total amount of SNPs for later normalization
        nbSNPs += p
    end
    # Normalization
    grm[:, :] /= nbSNPs

    close(outfile)
end
