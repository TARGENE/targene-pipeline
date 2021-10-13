using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Generate queries that can be used for TMLE, see "* 
                                     "https://github.com/olivierlabayle/GenesInteraction.jl/")

    @add_arg_table s begin
        "--mode", "-m"
            help = "Specifies how the queries are generated:"*
                   "The default,`frequency` will lead to a specification of the form "*
                   "MAJOR:MAJOR -> MAJOR:MINOR for each SNP"
            arg_type = String
            default = "frequency"
        "--out", "-o"
            help = "Out directory where the queries will be saved"
            arg_type = String
            default = "."
        "--sample-chrom-file", "-s"
            help = "Path to a sample chromosome (1-22) from the UKBB imputed directory. "*
                   "All chromosomes paths will be deduced from this file pattern."
            arg_type = String
            required = true
        "--threshold", "-t"
            arg_type = Float64
            help = "The threshold that will be used to hard call genotypes based on their probabilities"
            default = 0.9
        "--exclude", "-e"
            arg_type = String
            help = "The threshold that will be used to hard call genotypes based on their probabilities"
            required = false
        "asb-file"
            arg_type = String
            help = "Path to filtered allelic-specific binding SNPS"
            required = true
        "trans-actors-file"
            arg_type = String
            help = "Path to trans-acting SNPS"
            required = true
    end

    return parse_args(s)
end


parsed_args = parse_commandline()


generate_queries(parsed_args)
