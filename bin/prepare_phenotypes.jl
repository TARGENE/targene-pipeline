using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate queries that can be used for TMLE, see "* 
                                     "https://github.com/olivierlabayle/GenesInteraction.jl/")

    @add_arg_table s begin
        "binary-phenotypes"
            arg_type = String
            help = "Path to the GeneAtlas binary phenotypes file"
        "continuous-phenotypes"
            arg_type = String
            help = "Path to the GeneAtlas continuous phenotypes file"
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
            help = "Coma separated list"
            required = false
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


prepare_phenotypes(parsed_args)
