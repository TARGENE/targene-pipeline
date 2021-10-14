using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate queries that can be used for TMLE, see "* 
                                     "https://github.com/olivierlabayle/GenesInteraction.jl/")

    @add_arg_table s begin
        "input"
            arg_type = String
            help = "Path to the output of ukbconv containing data-fields"
            required = true
        "output"
            arg_type = String
            help = "Path where phenotypes will be saved with one column per phenotype."
            required = true
        "--save-cols", "-s"
            arg_type = String
            help = "Path where all extracted phenotypes names will be saved."
            required = false
        "--patterns", "-p"
            arg_type = String
            help = "A coma-separated list of patterns to match phenotypes with encoding "
                   "starting with any of the pattern."
            required = false
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


prepare_phenotypes(parsed_args)
