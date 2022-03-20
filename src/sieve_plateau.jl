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
The sample_grm comes from the gcta software. 
The process is as follows:
- Round between -1 and 1 as some elements may be beyond that limit
- Take 1 - this value to convert this quantity to a distance
- For each τ return if the distance between individuals is less or equal than τ
"""
function bit_distances(sample_grm, τs)
    distances = 1 .-  max.(min.(sample_grm, 1), -1)
    return convert(Matrix{Float32}, permutedims(distances) .<= τs)
end


default_τs(nτs;max_τ=2) = Float32[max_τ*(i-1)/(nτs-1) for i in 1:nτs]


function build_work_list(prefix, grm_ids; pval=0.05)
    dirname_, prefix_ = splitdir(prefix)
    dirname__ = dirname_ == "" ? "." : dirname_
    hdf5files = filter(
            x -> startswith(x, prefix_) && endswith(x, ".hdf5"), 
            readdir(dirname__)
    )
    hdf5files = [joinpath(dirname_, x) for x in hdf5files]

    influence_curves = Vector{Float32}[]
    n_obs = Int[]
    file_queryreport_pairs = Pair{String, Int}[]
    for hdf5file in hdf5files
        jldopen(hdf5file) do io
            qrs = haskey(io, "MACHINE") ? queryreports(io["MACHINE"]) : io["QUERYREPORTS"]
            for (qr_id, qr) in enumerate(qrs)
                if pvalue(ztest(qr)) <= pval
                    phenotype = string(qr.target_name)
                    sample_ids = parse.(Int, io["SAMPLE_IDS"][phenotype])
                    push!(influence_curves, align_ic(qr.influence_curve, sample_ids, grm_ids))
                    push!(n_obs, size(sample_ids, 1))
                    push!(file_queryreport_pairs, hdf5file => qr_id)
                end
            end
        end
    end
    return reduce(vcat, transpose(influence_curves)), n_obs, file_queryreport_pairs
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
    @views begin
        D_off_diag = transpose(influence_curves[:, 1:sample-1])
        D_diag = transpose(influence_curves[:, sample])
        return D_diag .* (2indicator[:, 1:sample-1] * D_off_diag .+ D_diag.* indicator[:, sample])
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
- variances: An Array of size (nτs, n_curves) where n_curves is the number of influence curves.
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
        variances .+= aggregate_variances(influence_curves, indicator, sample)
        start_idx = end_idx + 1
    end
    normalize!(variances, n_obs)
    return variances
end


function grm_rows_bounds(n_samples)
    bounds = Pair{Int, Int}[]
    start_idx = 1
    for sample in 1:n_samples
        # lower diagonal of the GRM are stored in a single vector 
        # that are accessed one row at a time
        end_idx = start_idx + sample - 1
        push!(bounds, start_idx => end_idx)
        start_idx = end_idx + 1
    end
    return bounds
end


function save_results(outfilename, τs, variances, std_errors, file_queryreport_pairs)
    jldopen(outfilename, "w") do io
        io["TAUS"] = τs
        io["VARIANCES"] = variances
        io["SOURCEFILE_REPORTID_PAIRS"] = file_queryreport_pairs
        io["STDERRORS"] = std_errors
    end
end


corrected_stderrors(variances, n_obs) =
    sqrt.(view(maximum(variances, dims=1), 1, :) ./ n_obs)


function sieve_variance_plateau(parsed_args)
    prefix = parsed_args["prefix"]
    pval = parsed_args["pval"]
    outfilename = parsed_args["out"]

    τs = default_τs(parsed_args["nb-estimators"];max_τ=parsed_args["max-tau"])
    grm, grm_ids = readGRM(parsed_args["grm-prefix"])

    influence_curves, n_obs, file_queryreport_pairs = build_work_list(prefix, grm_ids; pval=pval)
    variances = compute_variances(influence_curves, grm, τs, n_obs)
    std_errors = corrected_stderrors(variances, n_obs)

    save_results(outfilename, τs, variances, std_errors, file_queryreport_pairs)
end


