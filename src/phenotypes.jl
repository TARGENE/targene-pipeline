

function from_geneatlas(parsed_args)
    bridge = CSV.File(parsed_args["bridge"], header=["eid", "geneatlas_id"]) |> DataFrame
    # Loading in pure datafame objects seems to be fine with 16GB of memory
    binary_phenotypes  = CSV.File(parsed_args["binary-phenotypes"], drop=[:IID]) |> DataFrame
    continuous_phenotypes  = CSV.File(parsed_args["continuous-phenotypes"], drop=[:IID]) |> DataFrame
    merged = innerjoin(bridge, 
                        binary_phenotypes, 
                        on = :geneatlas_id => :FID)
    merged = innerjoin(merged,
                        continuous_phenotypes,
                        on=:geneatlas_id=> :FID)
    if haskey(parsed_args, "withdrawal-list")
        withdrawn = CSV.File(parsed_args["withdrawal-list"], header=["eid"]) |> DataFrame
        merged = filter(:eid => ∉(withdrawn.eid), merged)
    end

    CSV.write(parsed_args["output"], merged)
end

function prepare_phenotypes(parsed_args)
    if parsed_args["mode"] == "from-geneatlas"
        from_geneatlas(parsed_args)
    end
end