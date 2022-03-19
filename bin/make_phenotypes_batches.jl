using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate phenotypes batches to be passed to the TMLE step.")

    @add_arg_table s begin
        "phenotypes-file"
            arg_type = String
            help = "Path to a GeneAtlas phenotypes file"
        "--batch-size"
            arg_type = String
            help = "Number of phenotypes that will be contained in each batch passed for TMLE."*
                   "If 'max', a single batch is generated."
            default = "1"
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


tmle_phenotypes_batches(parsed_args)
