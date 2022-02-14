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


function update_variances!(variances, chunk_distance, inf_curve; row=1)
    grm_index = 1
    while grm_index < size(grm_chunk, 1)
        row += 1
        variances[k] += sum(skipmissing(chunk_distance[grm_index:grm_index+row-1] .* inf_curve[row] .* inf_curve[1:row]))
        grm_index += row
    end
    return (variances, row)
end


function grmparts_sizes(grmfiles)
    nparts = size(grmfiles, 1)
    grmpart_sizes = zeros(Int, nparts)
    for i in 1:nparts
        open(grmfiles[i]) do io
            grmpart_sizes[i] = read(io, Int)
        end
    end
    return grmpart_sizes
end

function bit_distance(grm_chunk_file, nτs)
    grmpart = grm_part(grm_chunk_file)
    distances = 1 .-  max.(min.(grmpart, 1), -1)
    τs = reshape([2/i for i in 1:nτs], 1, nτs)
    return distances .< τs
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
    influence_curves = zeros(Union{Float32, Missing},
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
            influence_curves[:, phenotype_id, query_id] = inf_curve
        end
    end
end


function sieve_plateau_variance(parsed_args)
    grm_ids = load_grm_ids(parsed_args["grm-ids"])
    results_file = jldopen(parsed_args["results"], "a+")

    nτs = parsed_args["nτs"]
    grmfiles = [joinpath(grmdir, file) for file in readdir(grmdir) if endswith(file, ".bin")]
    n_phenotypes = length(keys(results_file))

    variances = zeros(Float32, n_phenotypes, nτs)
    for grm_chunk_file in grmfiles
        chunk_distance = bit_distance(grm_chunk_file, nτs)
        for (phen_id, phenotype) in enumerate(keys(results_file))
            group = results_file[phenotype]
            sample_ids = group["sample_ids"]
            row = 1
            for queryreport in group["queryreports"]
                inf_curve = align_ic(queryreport.influence_curve, sample_ids, grm_ids)
                variances, row = update_variances!(variances, chunk_distance, inf_curve, τs; row=row)

            end
            group["variances"] = variances
        end
    end
end
