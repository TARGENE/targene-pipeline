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

function bit_distance(grm_chunk_file, nτs)
    grmpart = grm_part(grm_chunk_file)
    distances = 1 .-  max.(min.(grmpart, 1), -1)
    τs = reshape([2/i for i in 1:nτs], 1, nτs)
    return distances .<= τs
end


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
function compute_variances(influence_curves, nτs, grmfiles)
    n_phenotypes = size(influence_curves, 2)
    n_queries = size(influence_curves, 3)
    variances = zeros(Float32, n_phenotypes, n_queries, nτs)
    i_j_generator = generate_indices(grm_chunk_sizes(grmfiles), size(influence_curves, 1))
    first_i, first_j = 1, 1
    for grm_chunk_file in grmfiles
        last_i, last_j = take!(i_j_generator)
        chunk_distance = bit_distance(grm_chunk_file, nτs)
        for phenotype_idx in 1:n_phenotypes
            for query_idx in 1:n_queries
                inf_curve = influence_curves[:, phenotype_idx, query_idx]
                inf_pp = chunk_pairwise_inf_curve(inf_curve, first_i, first_j, last_i, last_j)
                variances[phenotype_idx, query_idx, :] .+= sum(chunk_distance .* inf_pp, dims=1)[1,:]
            end
        end
        if last_i == last_j
            first_i, first_j = last_i + 1, 1
        else
            first_i, first_j = last_i, last_j + 1
        end
    end
    return variances
end


pairwise_product(inf_curve, row; row_init=1, row_end=row) =
    inf_curve[row] .* inf_curve[row_init:row_end]

function chunk_pairwise_inf_curve(inf_curve, first_i, first_j, last_i, last_j)
    if first_j != 1
        pairwise_ic = pairwise_product(inf_curve, first_i, row_init=first_j)
        first_i += 1
        first_j = 1
    else
        pairwise_ic = Float32[]
    end

    for i in first_i:last_i-1
        pairwise_ic = vcat(pairwise_ic, pairwise_product(inf_curve, i))
    end

    (last_i == last_j) && return vcat(pairwise_ic, pairwise_product(inf_curve, last_i))

    return vcat(pairwise_ic, pairwise_product(inf_curve, last_i; row_end=last_j))
end


function generate_indices(sizes, n_samples)
    Channel{Tuple{Int, Int}}() do chnl
        sizes_ = copy(sizes)
        grm_chunk_size = popfirst!(sizes_)
        nb_elements = 0
        for i in 1:n_samples
            for j in 1:i
                nb_elements += 1
                if nb_elements == grm_chunk_size
                    put!(chnl, (i, j))
                    grm_chunk_size = popfirst!(sizes_)
                    nb_elements = 0
                end
            end
        end
    end
end