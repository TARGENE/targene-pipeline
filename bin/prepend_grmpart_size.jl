using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Prepends the size of a GRM part built by gcta64 "*
            "(https://yanglab.westlake.edu.cn/software/gcta/) to the bin file for faster reading.")

    @add_arg_table s begin
        "in-grmfile"
            help = "GRM part bin file as output by gcta64 --make-grm-parts"
            arg_type = String
            required = true
        "out-grmfile"
            help = "Outfile where the GRM part with size prepended information will be stored"
            arg_type = String
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


prepend_size(parsed_args)