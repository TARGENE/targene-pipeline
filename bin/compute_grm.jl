using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Computes the GRM given a list of bed files")

    @add_arg_table s begin
        "outdir"
            help = "A directory will be created and each line of the GRM will be represented as a file"
            required = true
        "bed-files"
            nargs = '+'
            arg_type = String
            help = "1 or more bed files to compute the GRM"
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


computeGRM(parsed_args)