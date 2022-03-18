

function from_geneatlas(parsed_args)
    bridge = CSV.File(parsed_args["bridge"], header=["SAMPLE_ID", "geneatlas_id"]) |> DataFrame
    # Loading in pure datafame objects seems to be fine with 16GB of memory
    phenotypes  = CSV.File(parsed_args["phenotypes-file"], drop=[:FID]) |> DataFrame

    merged = innerjoin(bridge, 
                       phenotypes, 
                       on = :geneatlas_id => :IID)
    
    select!(merged, Not(:geneatlas_id))

    if ~(parsed_args["phenotypes-list"] isa Nothing)
        filter_list = readlines(open(parsed_args["phenotypes-list"]))
        phen_list = filter(x-> x ∈ filter_list, names(phenotypes))
        merged = merged[!, vcat(["SAMPLE_ID"], phen_list)]
    end

    if ~(parsed_args["withdrawal-list"] isa Nothing)
        withdrawn = CSV.File(parsed_args["withdrawal-list"], header=["SAMPLE_ID"]) |> DataFrame
        merged = filter(:SAMPLE_ID => ∉(withdrawn.SAMPLE_ID), merged)
    end

    CSV.write(parsed_args["output"], merged)
end

function prepare_phenotypes(parsed_args)
    from_geneatlas(parsed_args)
end