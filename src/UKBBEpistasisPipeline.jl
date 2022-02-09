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
include("confounders.jl")
include("phenotypes.jl")
include("grm.jl")

# EXPORTS
#####################################################################

export generate_queries
export filter_asb
export filter_chromosome, merge_beds, adapt_flashpca
export prepare_phenotypes
export prepend_size

end
