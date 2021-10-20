# ADAPTING, MERGING AND FILTERING PRIOR TO PCA
#####################################################################


function merge_beds(parsed_args)
    merge_plink(parsed_args["input"]; des = parsed_args["output"])
end

issnp(x::String) = length(x) == 1

all_batches_ok(row, batchcols) = all(x == 1 for x in row[batchcols])

bounds(pos, lb, ub) =  (min(lb, pos - 10^4), max(ub, pos + 10^4))


function notin_ldblocks(row, ldblocks)
    if (chr=row.chromosome,) ∉ keys(ldblocks)
        return true
    else
        group = ldblocks[(chr=row.chromosome,)]
        inblock = any(b.lower_bound < row.position < b.upper_bound for b in eachrow(group))
        return ~inblock
    end
end


"""
    filter_chromosome(parsed_args)

The purpose of this method is to filter SNPS before applying a 
dimensionality reduction technique such as PCA that is later used
to extract population stratification compomemts.
We filter SNPs using quality control metrics from the following resource:
    - https://biobank.ndph.ox.ac.uk/showcase/refer.cgi?id=1955
"""
function filter_chromosome(parsed_args)

    qc_df = CSV.File(parsed_args["qcfile"]) |> DataFrame
    
    snp_data = SnpData(parsed_args["input"])
    # Load and redefine LD bounds
    ld_blocks = CSV.File(
        parsed_args["ld-blocks"];
        header=["rsid","chr","pos","LDblock_lower","LDblock_upper","LDblock_length","lower_bound","upper_bound"]) |> DataFrame
    ld_blocks.chr = string.(ld_blocks.chr)
    transform!(ld_blocks,
        :,
        [:pos, :lower_bound, :upper_bound] => ByRow(bounds) => [:lower_bound, :upper_bound],
        )
    ld_blocks = groupby(ld_blocks, :chr)

    # Remove SNP's with MAF < maf-threshold
    maf_threshold = parsed_args["maf-threshold"]
    snp_data.snp_info[!, "MAF"] = SnpArrays.maf(snp_data.snparray)
    mafpassed = filter(:MAF => >=(maf_threshold), snp_data.snp_info)
    
    # Remove LD regions specified by ld_blocks
    ld_pruned = filter(x -> notin_ldblocks(x, ld_blocks), mafpassed)

    # The QC file contains information on fully genotyped SNPS
    # We only keep those
    fully_genotyped_snps = innerjoin(
        ld_pruned, 
        qc_df, 
        on = :snpid => :rs_id,
        makeunique = true
        )

    # If an RSID appears multiple times, it is because it has 
    # more than 2 possible alleles: we remove them 
    # (why? maybe because the PCA then cannot tackle them)
    duplicate_rsids = Set(fully_genotyped_snps.snpid[nonunique(fully_genotyped_snps, ["snpid"])])
    biallelic = filter(:snpid=>∉(duplicate_rsids), fully_genotyped_snps)

    # Keep only actual SNPs and not other kinds of variants
    actual_snps = subset(biallelic, :allele1 => ByRow(issnp), :allele2 => ByRow(issnp))

    # All batches pass QC 
    batch_cols = [x for x in names(actual_snps) if occursin("Batch", x)]
    batches_ok = filter(row -> all_batches_ok(row, batch_cols), actual_snps)

    # Assayed in both genotyping arrays
    final = filter(:array => ==(2), batches_ok)
    
    rsids = Set(final.snpid)
    SnpArrays.filter(parsed_args["input"]; des=parsed_args["output"], f_snp = x -> x[:snpid] ∈ rsids)
end


function adapt_flashpca(parsed_args)
    # I think FID and IID are just duplicates
    pcs = CSV.File(parsed_args["input"], drop=["IID"]) |> DataFrame
    rename!(pcs, :FID => :SAMPLE_ID)
    CSV.write(parsed_args["output"], pcs)
end