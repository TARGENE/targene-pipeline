GRMIDs(file) = CSV.File(file, 
                        header=["FAMILY_ID", "SAMPLE_ID"], 
                        select=["SAMPLE_ID"]) |> DataFrame

function readGRM(prefix)
    ids = GRMIDs(string(prefix, ".id"))
    n = size(ids, 1)
    grm_size = n*(n + 1) ÷ 2
    GRM = mmap(string(prefix, ".bin"), Vector{Float32}, grm_size)

    return GRM, ids
end