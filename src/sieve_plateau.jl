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

function bit_distances(sample_grm, nτs)
    distances = 1 .-  max.(min.(sample_grm, 1), -1)
    τs = reshape(default_τs(nτs), 1, nτs)
    return distances .<= τs
end

default_τs(nτs) = [2/i for i in 1:nτs]

"""
    build_influence_curves(results_file, grm_ids)

Load and build an influence curves 3-dimensional array with sizes:
    - N samples 
    - N phenotypes 
    - N queries

The ordering of the samples in the built influence_curves 
matrix will match the GRM sample_ids.
"""
function build_influence_curves(results_file, grm_ids)
    phenotypes = keys(results_file)
    any_phenotype = first(phenotypes)
    any_queryreports = results_file[any_phenotype]["queryreports"]
    influence_curves = zeros(Float32,
        size(grm_ids, 1),
        size(phenotypes, 1),
        size(any_queryreports, 1)
        )
    fill_influence_curves!(influence_curves, results_file, phenotypes, grm_ids)
    return phenotypes, influence_curves
end


function fill_influence_curves!(influence_curves, results_file, phenotypes, grm_ids)
    for (phenotype_id, phenotype) in enumerate(phenotypes)
        phenotype_group = results_file[phenotype]
        sample_ids = phenotype_group["sample_ids"]
        queryreports = results_file[phenotype]["queryreports"]
        for (query_id, queryreport) in enumerate(queryreports)
            inf_curve = align_ic(queryreport.influence_curve, sample_ids, grm_ids)
            influence_curves[:, phenotype_id, query_id] = coalesce.(inf_curve, 0)
        end
    end
end


"""
    compute_variances(influence_curves, nτs, grm_files)

Each influence_curve vector is associated with a phenotype query.
This function compute the variance estimates at each τ for all queries.
"""
function compute_variances(influence_curves, grm, nτs)
    n_samples, n_phenotypes, n_queries = size(influence_curves)
    variances = zeros(Float32, n_phenotypes, n_queries, nτs)
    for sample_idx in 1:n_samples
        bit_distances_ = bit_distances(grm[:, sample_idx], nτs)
        for phenotype_idx in 1:n_phenotypes
            for query_idx in 1:n_queries
                inf_curve = influence_curves[:, phenotype_idx, query_idx]
                variances[phenotype_idx, query_idx, :] .+= sum(inf_curve[sample_idx] .* bit_distances_ .* inf_curve, dims=1)[1, :]
            end
        end
    end
    return variances
end

function update_results_file!(results_file, phenotypes, variances)
    for (phenotype_idx, phenotype) in enumerate(phenotypes)
        results_file[phenotype]["variances"] = variances[phenotype_idx, :, :]
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

    nτs = parsed_args["nb-estimators"]
    grm, grm_ids = readGRM(parsed_args["grm-prefix"])
    phenotypes, inf_curves = build_influence_curves(results_file, grm_ids)
    variances = compute_variances(inf_curves, grm, nτs)
    #monotone_variances, radial_variances = make_monotone_and_smooth(variances)

    update_results_file!(results_file, phenotypes, variances)

    close(results_file)
end