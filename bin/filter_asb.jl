using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Merge and filter allelic-specific binding SNPs output by "* 
                                     "https://git.ecdf.ed.ac.uk/tfomics/binding-variants/baal-process-vdr")

    @add_arg_table s begin
        "--mode", "-m"
            help = "The rrule used to filter the SNPs: \n"*
                   " - default: `isASB=true` are kept"
            arg_type = String
            default = "default"
        "--out", "-o"
            help = "Where to save thre filtered ASBs"
            required = true
        "asb-files"
            nargs = '+'
            arg_type = String
            help = "1 or more files that will be merged and filtered"
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


filter_asb(parsed_args)



