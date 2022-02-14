
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


load_grm_ids(grm_id_file) = 
    CSV.File(grm_id_file,
            select=[2],
            header=["FAMILY_ID", "SAMPLE_ID"]) |> DataFrame