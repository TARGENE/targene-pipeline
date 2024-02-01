include { longest_prefix } from './utils.nf'

process GRMPart {
    container "olivierlabayle/tl-core:cvtmle"
    label "bigmem"
    label "multithreaded"

    input:
        path bedfiles
        val nparts
        val part_id

    output:
        path "GRM*.grm.*"
    
    script:
        base = bedfiles.first().getName().split("\\.")[0]
        "gcta64 --bfile $base --make-grm-part $nparts $part_id --thread-num ${task.cpus} --out GRM"
}

process AggregateGRM {
    publishDir "$params.OUTDIR/GRM", mode: 'symlink'

    input:
        path grm_files

    output:
        path "GRM.id", emit: grm_ids
        path "GRM.bin", emit: grm_matrix

    script:
        """
        cat *.grm.id > GRM.id
        cat *.grm.bin > GRM.bin
        """
}

process SVP {
    container "olivierlabayle/targeted-estimation:argparse"
    publishDir "$params.OUTDIR", mode: 'symlink'

    input:
        path hdf5_results
        path GRM_ids
        path GRM_matrix
        val n_estimators
        val max_tau
        val estimator_key
        val verbosity

    output:
        path "svp.hdf5"
    
    script:
        hdf5_prefix = longest_prefix(hdf5_results)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --sysimage=/TargetedEstimation.jl/TMLESysimage.so --project=/TargetedEstimation.jl --startup-file=no /TargetedEstimation.jl/tmle.jl svp \
        $hdf5_prefix \
        --n-estimators=$n_estimators \
        --max-tau=$max_tau \
        --estimator-key=$estimator_key \
        --verbosity=$verbosity
        """
}