module UKBBEpistasisPipeline

using MKL
using DataFrames
using CSV
using TOML
using BGEN
using SnpArrays
using Mmap
using JLD2
using TMLE

# INCLUDES
#####################################################################

include("asb_snps.jl")
include("queries_generation.jl")
include("confounders.jl")
include("phenotypes.jl")
include("grm.jl")
include("sieve_plateau.jl")

# EXPORTS
#####################################################################

export generate_queries
export filter_asb
export filter_chromosome, merge_beds, adapt_flashpca
export prepare_phenotypes
export grm_from_gcta
export sieve_variance_plateau

end
