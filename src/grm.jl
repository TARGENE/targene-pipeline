
function fillgrmpart!(grmpart, io)
    while !eof(io)
        push!(grmpart, read(io, Float32))
    end
end

function ReagGRMPart(file)
    io = open(file)
    grmpart = Float32[]
    fillgrmpart!(grmpart, io)
    return DataFrame(RELATIONSHIPS=grmpart)
end


function grm_parts_to_arrow(parsed_args)
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
                part_ids = CSV.File(joinpath(dir, file), header=["FAMILY_ID", "ID"]) |> DataFrame
                ids = vcat(ids, part_ids)
            end
        end
    end
    # Write arrow GRM
    grm_parts = Tables.partitioner(ReagGRMPart, binfiles)
    Arrow.write(string(outprefix, ".arrow"), grm_parts)
    # Write ids
    CSV.write(string(outprefix, ".ids.csv"), ids)
end

