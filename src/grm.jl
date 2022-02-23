
function fillgrmpart!(grmpart, io)
    while !eof(io)
        push!(grmpart, read(io, Float32))
    end
end

function grm_part_from_gcta(file)
    io = open(file)
    grmpart = Float32[]
    fillgrmpart!(grmpart, io)
    return grmpart
end


"""
    prepend_size(parsed_args)

This unfortunate procedure is intended for reading speed. If the size of the 
array on disk is known in advance then reading is improved by 100 folds.
"""
function prepend_size(parsed_args)
    grmpart = grm_part_from_gcta(parsed_args["in-grmfile"])
    open(parsed_args["out-grmfile"], "w") do io
        write(io, size(grmpart, 1))
        write(io, grmpart)
    end
end


function grm_part(grmpart_file)
    io = open(grmpart_file)
    grmpart_size = read(io, Int)
    GRMPart = Vector{Float32}(undef, grmpart_size)
    read!(io, GRMPart)
    return GRMPart
end


function grm_chunk_sizes(grm_files)
    chunk_sizes = zeros(Int, size(grm_files, 1))
    for (i, grm_chunk_file) in enumerate(grm_files)
        open(grm_chunk_file) do io
            chunk_sizes[i] = read(io, Int)
        end
    end
    return chunk_sizes
end

load_grm_ids(grm_id_file) = 
    CSV.File(grm_id_file,
            select=[2]) |> DataFrame


function update_grm!(GRM, grm_chunk, index, i, binfiles)
    if i + index - 1 <= length(grm_chunk)
        row_chunk = grm_chunk[index:index + i - 1]
        index += i
    else
        row_chunk = grm_chunk[index:end]
        index = i - length(row_chunk) + 1
        @info string("Loading file: ", binfiles[1])
        grm_chunk = grm_part_from_gcta(popfirst!(binfiles))
        index - 1 < length(grm_chunk) || throw(BoundsError("A grm chunk has fewer than expected number of elements"))
        row_chunk = vcat(row_chunk, grm_chunk[1:index-1])
    end
    GRM[1:i, i] = row_chunk
    return grm_chunk, index
end

function grm_from_gcta(parsed_args)
    inprefix = parsed_args["inprefix"]
    outprefix = parsed_args["outprefix"]
    dir, baseprefix = splitdir(inprefix)

    binfiles = String[]
    ids = DataFrame()

    dir_ = dir == "" ? "." : dir
    # The sorting ensures that grm parts are ordered
    for file in sort(readdir(dir_))
        # Read files only pertaining to the GRM parts
        if contains(file, ".part_")
            if endswith(file, "grm.bin") && contains(file, ".part_")
                push!(binfiles, joinpath(dir, file))
            elseif endswith(file, "grm.id") && startswith(file, baseprefix)
                part_ids = CSV.File(joinpath(dir, file), header=["FAMILY_ID", "SAMPLE_ID"]) |> DataFrame
                ids = vcat(ids, part_ids)
            end
        end
    end
    # Write unified GRM
    n = size(ids, 1)
    index = 1
    @info string("Loading file: ", binfiles[1])
    grm_chunk = grm_part_from_gcta(popfirst!(binfiles))
    GRM = mmap(string(outprefix, ".bin"), Matrix{Float32}, (n, n))
    @inbounds for i in 1:n
        grm_chunk, index = update_grm!(GRM, grm_chunk, index, i, binfiles)
    end
    # Write ids
    CSV.write(string(outprefix, ".ids.csv"), ids)

end

function readGRM(prefix)
    ids = load_grm_ids(string(prefix, ".ids.csv"))
    n = size(ids, 1)
    GRM = mmap(string(prefix, ".bin"), Matrix{Float32}, (n, n))
    return Symmetric(GRM), ids
end