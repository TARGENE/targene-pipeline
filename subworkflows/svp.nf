include { SVP; AggregateGRM; GRMPart } from '../modules/svp.nf'

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