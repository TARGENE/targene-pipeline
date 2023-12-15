process GRMPart {
    container "olivierlabayle/tl-core:0.6"
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
    container "olivierlabayle/targeted-estimation:cv_tmle"
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
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargetedEstimation.jl --startup-file=no /opt/bin/tmle sieve-variance-plateau \
        $hdf5_prefix \
        --n-estimators=$n_estimators \
        --max-tau=$max_tau \
        --estimator-key=$estimator_key \
        --verbosity=$verbosity
        """
}


workflow SVPWorkflow {
    take:
        hdf5_result
        iid_genotypes
        n_svp_estimators
        max_svp_threshold
        svp_estimator_key
        grm_n_splits
        verbosity
    
    main:
        grm_parts = Channel.from( 1..grm_n_splits )
        GRMPart(iid_genotypes.collect(), grm_n_splits, grm_parts)
        AggregateGRM(GRMPart.out.collect())
        // Sieve estimation
        SVP(
            hdf5_result.collect(), 
            AggregateGRM.out.grm_ids, 
            AggregateGRM.out.grm_matrix,
            n_svp_estimators,
            max_svp_threshold,
            svp_estimator_key,
            verbosity
        )
}