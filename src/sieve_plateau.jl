function align_ic(ic, sample_ids, grm_ids)
    leftjoin!(
        grm_ids, 
        DataFrame(IC=ic, SAMPLE_ID=sample_ids),
        on=:SAMPLE_ID
    )
    aligned_ic = grm_ids.IC
    select!(grm_ids, Not(:IC))
    return coalesce.(aligned_ic, 0)
end

"""
    bit_distances(sample_grm, nτs)

Returns a matrix of shape (n_samples,nτs) where n_samples is the 
size of sample_grm.
The sample_grm comes from the gcta software which is in fact GRM/2. 
The process is thus as follows:
- Multiply by 2 to get back to the GRM
- Round between -1 and 1 as some elements may be beyond that limit
- Take 1 - this value to covnert this quantity to a distance
- For each τ return if the distance between individuals is less or equal than τ
"""
function bit_distances(sample_grm, τs)
    distances = 1 .-  max.(min.(2sample_grm, 1), -1)
    return convert(Matrix{Float32}, permutedims(distances) .<= τs)
end

default_τs(nτs) = [2/i for i in 1:nτs]

function build_work_list(results_file, grm_ids; pval=0.05)
    influence_curves = Vector{Float32}[]
    n_obs = Int[]
    phenotype_query_pairs = Pair{String, Int}[]
    for phenotype in keys(results_file)
        for (query_idx, queryreport) in enumerate(results_file[phenotype]["queryreports"])
            if pvalue(ztest(queryreport)) <= pval
                sample_ids = parse.(Int, results_file[phenotype]["sample_ids"])
                push!(influence_curves, align_ic(queryreport.influence_curve, sample_ids, grm_ids))
                push!(n_obs, size(sample_ids, 1))
                push!(phenotype_query_pairs, phenotype=>query_idx)
            end
        end
    end
    return reduce(vcat, transpose(influence_curves)), n_obs, phenotype_query_pairs
end



"""
    normalize(variances, n_observations)

Divides the variance estimates by the effective number of observations 
used for each phenotype at estimation time.
"""
normalize!(variances, n_observations) = 
    variances ./= permutedims(n_observations)


"""
    aggregate_variances(influence_curves, indicator, sample)

This function computes the sum for a single index i, see also `compute_variances`.
As the GRM is symetric it is performed as : 
    2 times off-diagonal elements with j < i + diagonal term 
and this for all τs.
"""
function aggregate_variances(influence_curves, indicator, sample)
    D_off_diag = transpose(influence_curves[:, 1:sample-1])
    D_diag = transpose(influence_curves[:, sample])
    return D_diag .* (2indicator[:, 1:sample-1] * D_off_diag .+ D_diag.* indicator[:, sample])
end

"""
    compute_variances(influence_curves, nτs, grm_files)

An overall variance estimate for a distance function d, a threshold τ 
and influence curve D is given by:
            σ̂ = 1/n ∑ᵢ∑ⱼ1(dᵢⱼ ≤ τ)DᵢDⱼ
              = 1/n 2∑ᵢDᵢ∑ⱼ<ᵢ1(dᵢⱼ ≤ τ)DᵢDⱼ + 1(dᵢᵢ ≤ τ)DᵢDᵢ

This function computes those variance estimates at each τ for all phenotypes
and queries.

# Arguments:
- influence_curves: Array of size (n_samples, n_queries, n_phenotypes)
- grm: Vector containing the lower elements of the GRM ie: n(n+1)/2 elements
- τs: 1 row matrix containing the distance thresholds between individuals
- n_obs: Vector of size (n_phenotypes,), containing the number of effective 
observations used during estimation

# Returns:
- variances: An Array of size (nτs, n_queries, n_phenotypes)
"""
function compute_variances(influence_curves, grm, τs, n_obs)
    n_curves, n_samples = size(influence_curves)
    variances = zeros(Float32, length(τs), n_curves)
    start_idx = 1
    for sample in 1:n_samples
        # lower diagonal of the GRM are stored in a single vector 
        # that are accessed one row at a time
        end_idx = start_idx + sample - 1
        sample_grm = view(grm, start_idx:end_idx)
        indicator = bit_distances(sample_grm, τs)
        res = aggregate_variances(influence_curves, indicator, sample)
        variances .+= res
        start_idx = end_idx + 1
    end
    normalize!(variances, n_obs)
    return variances
end

function update_results_file!(results_file, variances, phenotype_query_pairs)
    for (curve_id, (phenotype, query_idx)) in enumerate(phenotype_query_pairs)
        if !haskey(results_file[phenotype], "sieve_variances")
            sieve_group = JLD2.Group(results_file[phenotype], "sieve_variances")
        else
            sieve_group = results_file[phenotype]["sieve_variances"]
        end
        sieve_group[string(query_idx)] = variances[:, curve_id]
    end
end

function sieve_variance_plateau(parsed_args)
    results_file = jldopen(parsed_args["results"], "a+")

    τs = default_τs(parsed_args["nb-estimators"])
    grm, grm_ids = readGRM(parsed_args["grm-prefix"])
    influence_curves, n_obs, phenotype_query_pairs = build_work_list(results_file, grm_ids; pval=0.05)
    variances = compute_variances(influence_curves, grm, τs, n_obs)
    #monotone_variances, radial_variances = make_monotone_and_smooth(variances)

    update_results_file!(results_file, variances, phenotype_query_pairs)

    close(results_file)
end

"""
make_monotone_and_smooth(variances)

For each phenotype/query, build monotonically increasing variances estimates and 
smooth them with radial basis functions.

# Arguments:
- variances: an Array with shape: (n_phenotypes, n_queries, nτ)

"""
function make_monotone_and_smooth(variances)
    n_phenotypes, n_queries, nτ = size(variances)
    monotone_variances = zeros(n_phenotypes, n_queries, nτ)
    radial_variances = zeros(n_phenotypes, n_queries, nτ)
    for query_idx in 1:n_queries
        for phenotype_idx in 1:n_phenotypes
            monotone_variances[phenotype_idx, query_idx, :] = 
                pooled_pava_isotonic_regression(variances[phenotype_idx, query_idx, :])
            radial_variances[phenotype_idx, query_idx, :] =
                radial_basis_interpolation(monotone_variances[phenotype_idx, query_idx, :])
        end
    end
    monotone_variances, radial_variances
end


