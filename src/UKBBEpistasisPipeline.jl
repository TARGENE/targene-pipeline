module UKBBEpistasisPipeline

using DataFrames
using CSV
using TOML
using BGEN
using SnpArrays

# INCLUDES
#####################################################################

include("asb_snps.jl")
include("queries_generation.jl")
include("confounders_preparation.jl")
include("phenotypes.jl")


# EXPORTS
#####################################################################

export generate_queries
export filter_asb
export filter_chromosome, merge_beds

end
