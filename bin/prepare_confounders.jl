using UKBBEpistasisPipeline
using ArgParse


function parse_commandline()
    s = ArgParseSettings(description="Utilities to prepare genotypes for later PCA extraction.")

    @add_arg_table s begin
        "--input", "-i"
            help = "If merge is called, a common prefix to the bed files should be provided."*
                "If filter is called, a path to the chromosome .bed file should be provided"
            arg_type = String
            required = true
        "--output", "-o"
            help = "Path to the output file, if you use merge after filter you may want to use a common prefix"
            arg_type = String
            required = true
        "--qcfile"
            help = "Path to the UKBiobank ukb_snp_qc.txt"
            arg_type = String
        "--maf-threshold", "-t"
            help = "SNPs with MAF lower than this value will be filtered out"
            arg_type = Float64
            default = 0.01
            required = false
        "--ld-blocks", "-l"
            help = "Path to the LD block definition around the SNPs of interest. "*
                   "Issued from an LD analysis."
            arg_type = String
            required = false
        "filter"
            help = "Filter the SNPS from a PLINK .bed file based on quality metrics"
            action = :command
        "merge"
            help = "Merge multiple PLINK .bed files together"
            action = :command
        "adapt"
            help = "Adapts the csv output of flashpca to what is expected by the TMLE step"
            action = :command

    end

    return parse_args(s)
end


parsed_args = parse_commandline()

if haskey(parsed_args, "merge") == true
    merge_beds(parsed_args)
elseif haskey(parsed_args, "filter") == true
    filter_chromosome(parsed_args)
elseif haskey(parsed_args, "adapt")
    adapt_flashpca(parsed_args)
end
