function align_ic(ic, sample_ids, grm_ids)
    leftjoin!(
        grm_ids, 
        DataFrame(IC=ic, SAMPLE_ID=sample_ids),
        on=:SAMPLE_ID
    )
    aligned_ic = grm_ids.IC
    select!(grm_ids, Not(:IC))
    return aligned_ic
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
    return convert(Matrix{Float32}, distances .<= τs)
end

default_τs(nτs) = reshape([2/i for i in 1:nτs], 1, nτs)

"""
    build_influence_curves(results_file, grm_ids)

Load and build an influence curves 3-dimensional array with sizes:
    - N samples 
    - N queries
    - N phenotypes 
    
The ordering of the samples in the built influence_curves 
matrix will match the GRM sample_ids.
"""
function build_influence_curves(results_file, grm_ids)
    phenotypes = keys(results_file)
    any_phenotype = first(phenotypes)
    any_queryreports = results_file[any_phenotype]["queryreports"]
    n_observations = zeros(Int, length(phenotypes))
    influence_curves = zeros(Float32,
        size(grm_ids, 1),
        size(any_queryreports, 1),
        size(phenotypes, 1),
        )
    fill_influence_curves!(influence_curves, n_observations, results_file, phenotypes, grm_ids)
    return phenotypes, influence_curves, n_observations
end


function fill_influence_curves!(influence_curves, n_observations, results_file, phenotypes, grm_ids)
    for (phenotype_id, phenotype) in enumerate(phenotypes)
        phenotype_group = results_file[phenotype]
        sample_ids = phenotype_group["sample_ids"]
        n_observations[phenotype_id] = size(sample_ids, 1)
        queryreports = results_file[phenotype]["queryreports"]
        for (query_id, queryreport) in enumerate(queryreports)
            inf_curve = align_ic(queryreport.influence_curve, sample_ids, grm_ids)
            influence_curves[:, query_id, phenotype_id] = coalesce.(inf_curve, 0)
        end
    end
end

"""
    aggregate_variance(D, indicator, i)

This function computes the sum for a single index i, see also `compute_variances`.
As the GRM is symetric it is performed as : 
    2 times off-diagonal elements with j < i + diagonal term 
and this for all τs.
"""
aggregate_variance(D, indicator) =
    @views D[end].*(2sum(indicator[1:end-1, :] .* D[1:end-1], dims=1)[1, :] .+ D[end] .* indicator[end, :])

"""
    normalize(variances, n_observations)

Divides the variance estimates by the effective number of observations 
used for each phenotype at estimation time.
"""
function normalize!(variances, n_observations)
    n_phenotypes = size(variances, 3)
    for i in 1:n_phenotypes
        variances[:, :, i] ./= n_observations[i]
    end
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
    n_samples, n_queries, n_phenotypes = size(influence_curves)
    variances = zeros(Float32, length(τs), n_queries, n_phenotypes)
    start_idx = 1
    @inbounds for sample in 1:n_samples
        # lower diagonal of the GRM are stored in a single vector 
        # that are accessed one row at a time
        end_idx = start_idx + sample - 1
        sample_grm = view(grm, start_idx:end_idx)
        indicator = UKBBEpistasisPipeline.bit_distances(sample_grm, τs)
        @sync for query_idx in 1:n_queries
            for phenotype_idx in 1:n_phenotypes
                @spawn variances[:, query_idx, phenotype_idx] .+= 
                    UKBBEpistasisPipeline.aggregate_variance(view(influence_curves, 1:sample, query_idx, phenotype_idx), indicator)
            end
        end
        start_idx = end_idx + 1
    end
    normalize!(variances, n_obs)
    return variances
end

function update_results_file!(results_file, phenotypes, variances)
    for (phenotype_idx, phenotype) in enumerate(phenotypes)
        results_file[phenotype]["variances"] = variances[:, :, phenotype_idx]
    end
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


function sieve_variance_plateau(parsed_args)
    results_file = jldopen(parsed_args["results"], "a+")

    τs = default_τs(parsed_args["nb-estimators"])
    grm, grm_ids = readGRM(parsed_args["grm-prefix"])
    phenotypes, inf_curves, n_obs = build_influence_curves(results_file, grm_ids)
    variances = compute_variances(inf_curves, grm, τs, n_obs)
    #monotone_variances, radial_variances = make_monotone_and_smooth(variances)

    update_results_file!(results_file, phenotypes, variances)

    close(results_file)
end

