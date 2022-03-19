

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

function write_phenotype_batch(phenotypes, batch, batch_start, batch_end)
    batch_file = batch_filename(batch)
    open(batch_file, "w") do io
        for i in batch_start:batch_end
            println(io, phenotypes[i])
        end
    end
end

batch_filename(batch) = string("phenotypes_batch_", batch, ".csv")

function tmle_phenotypes_batches(parsed_args)
    phenotypes = filter(
        !=("SAMPLE_ID"), 
        split(readline(open(parsed_args["phenotypes-file"])), ",")
    )
    n_phenotypes = size(phenotypes, 1)
    batch_size = parsed_args["batch-size"]
    batch_size = batch_size == "max" ? size(phenotypes, 1) : parse(Int, batch_size)
    
    batch = 1
    batch_start = 1
    batch_end = batch_start + batch_size - 1
    while batch_end <= n_phenotypes
        write_phenotype_batch(phenotypes, batch, batch_start, batch_end)
        batch += 1
        batch_start = batch_end + 1
        batch_end = batch_start + batch_size - 1
    end

    if batch_start <= n_phenotypes
        write_phenotype_batch(phenotypes, batch, batch_start, n_phenotypes)
    end

end