using UKBBEpistasisPipeline
using ArgParse


using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Aggregate parts of the GRM built by gcta64 "*
            "(https://yanglab.westlake.edu.cn/software/gcta/) to Apache Arrow format.")

    @add_arg_table s begin
        "inprefix"
            help = "Prefix to the GRM parts output by gcta64 --make-grm-parts"
            arg_type = String
            required = true
        "outprefix"
            help = "Output prefix for both Arrow aggregated GRM and associated sample ids."
            arg_type = String
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


grm_from_gcta(parsed_args)