include { longest_prefix } from './utils.nf'

process GRMPart {
    label "bigmem"
    label "multithreaded"
    label 'targenecore_image'

    input:
        path bedfiles
        val nparts
        val part_id

    output:
        path "GRM*.grm.*"
    
    script:
        // Remove the last "." from the prefix otherwise gcta will not find the files
        input_prefix = longest_prefix(bedfiles)[0..-2]
        "gcta64 --bfile ${input_prefix} --make-grm-part ${nparts} ${part_id} --thread-num ${task.cpus} --out GRM"
}

process AggregateGRM {
    publishDir "${params.OUTDIR}/GRM", mode: 'symlink'

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
    label 'targenecore_image'
    publishDir "${params.OUTDIR}", mode: 'symlink'

    input:
        path hdf5_results
        path GRM_ids
        path GRM_matrix

    output:
        path "svp.hdf5"
    
    script:
        hdf5_prefix = longest_prefix(hdf5_results)
        """
        TEMPD=\$(mktemp -d)
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TargeneCore.jl --startup-file=no --sysimage=/TargeneCore.jl/TargeneCoreSysimage.so /TargeneCore.jl/targenecore.jl svp \
        ${hdf5_prefix} \
        --n-estimators=${params.NB_SVP_ESTIMATORS} \
        --max-tau=${params.MAX_SVP_THRESHOLD} \
        --estimator-key=${params.ESTIMATOR_KEY} \
        --verbosity=${params.VERBOSITY}
        """
}
