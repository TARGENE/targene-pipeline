using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate queries that can be used for TMLE, see "* 
                                     "https://github.com/olivierlabayle/GenesInteraction.jl/")

    @add_arg_table s begin
        "phenotypes-file"
            arg_type = String
            help = "Path to a GeneAtlas phenotypes file"
        "bridge"
            arg_type = String
            help = "Path to the bridge file between our UKBB project and the GeneAtlas."
        "output"
            arg_type = String
            help = "Path where all phenotypes will be saved"
        "--withdrawal-list"
            arg_type = String
            help = "Path to participants withdrawal list"
            required = false
        "--phenotypes-list"
            arg_type = String
            help = "A path to a file containing a list of phenotypes of interest. One per line."
            required = false
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


prepare_phenotypes(parsed_args)
