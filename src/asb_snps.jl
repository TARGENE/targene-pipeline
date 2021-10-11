# FILTERING ASB SNPS
#####################################################################

function default_filter_asb(asb_files, out)
    filtered_asbs = DataFrame()
    for file in asb_files
        asb_df = CSV.File(file) |> DataFrame
        filtered_asbs = vcat(filtered_asbs, filter(x -> x.isASB === true, asb_df))
    end

    CSV.write(out, filtered_asbs)
end


function filter_asb(parsed_args)
    if parsed_args["mode"] == "default"
        default_filter_asb(parsed_args["asb-files"], parsed_args["out"])
    else
        throw(ArgumentError("This method is not implemented."))
    end
end