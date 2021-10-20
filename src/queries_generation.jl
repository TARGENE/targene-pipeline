# GENERATING QUERIES
#####################################################################

function get_minor_major(parsed_args)
    # Load SNPS
    trans_actors = CSV.File(parsed_args["trans-actors-file"], select=[:ID,:CHROM]) |> DataFrame
    asbs = CSV.File(parsed_args["asb-file"], select=[:ID,:CHROM]) |> DataFrame
    
    # Get base file name
    basedir, chrfile = splitdir(parsed_args["sample-chrom-file"])
    re = r"chr[1-9]+"
    
    extra_info = DataFrame(ID = String[], MAJOR = String[], MINOR = String[], CHROMPATH = String[])
    
    allsnps = vcat(trans_actors, asbs)
    # Remove potentially unwanted snps
    if ~(parsed_args["exclude"] isa Nothing)
        snps_to_remove = collect(CSV.File(parsed_args["exclude"], header=["SNPS_TO_REMOVE"]).SNPS_TO_REMOVE)
        filter!(:ID => x -> x ∉ snps_to_remove, allsnps)
    end

    # Load by chromosome file
    chrom_groups = groupby(allsnps, :CHROM)
    for (chr, group) in pairs(chrom_groups)
        chrompath = joinpath(basedir, replace(chrfile, re => chr.CHROM))
        b = BGEN.Bgen(chrompath)
        for (rsid, _) in eachrow(group)

            v = variant_by_rsid(b, rsid)
            minor_allele_dosage!(b, v)

            minor = minor_allele(v)
            major = major_allele(v)

            push!(extra_info, (rsid, major, minor, chrompath))
        end
        
    end

    return innerjoin(trans_actors, extra_info, on=:ID), innerjoin(asbs, extra_info, on=:ID)
end

function major_minor_queries(parsed_args)
    threshold = parsed_args["threshold"]
    # Load SNPS
    trans_actors, asbs = get_minor_major(parsed_args)
    # Loop over cartesian product of asb and trans
    merged = crossjoin(asbs, trans_actors, makeunique=true)
    for row in eachrow(merged)
        query = Dict(
            "threshold" => threshold,
            "SNPS" => Dict(
                row.ID => row.CHROMPATH,
                row.ID_1 => row.CHROMPATH_1
            ),
            "QUERY" => Dict(
                row.ID => row.MAJOR*row.MAJOR*" -> "*row.MAJOR*row.MINOR,
                row.ID_1 => row.MAJOR_1*row.MAJOR_1*" -> "*row.MAJOR_1*row.MINOR_1
            )
        )

        outfile = joinpath(parsed_args["out"], "query_$(row.ID)_$(row.ID_1).toml")
        open(outfile, "w") do io
            TOML.print(io, query)
        end
    end

end


function generate_queries(parsed_args)
    if parsed_args["mode"] == "frequency"
        major_minor_queries(parsed_args)
    else
        throw(ArgumentError("Not implemented yet."))
    end
end