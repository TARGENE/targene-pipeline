
function columnnames(filename::String, patterns)
    columns = []
    for row in CSV.File(filename)
        
    end
end

function prepare_phenotypes(parsed_args)
    patterns = parsed_args["subset"]
    
end