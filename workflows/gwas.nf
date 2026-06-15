include { LocoPCA } from './pca.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'
include { IIDGenotypes } from '../subworkflows/confounders.nf'
include { mergeBEDS } from '../modules/confounders.nf'
include { ReportWorkflow } from '../subworkflows/report.nf'
include { subsetBED; denseBED } from '../modules/extract_variants.nf'

workflow GWAS {
    // Define Parameters
    bed_files = channel.fromFilePairs("$params.BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
    estimands_file = channel.value(file("$params.ESTIMANDS_CONFIG"))
    estimator_config = channel.fromPath(EstimatorsConfig.create(params.ESTIMATORS_CONFIG, params.OUTDIR))

    // Loco PCA
    LocoPCA()
    
    // Estimation Inputs
    if (params.SUBSET_FILE != "NO_SUBSET_FILE") {
        subset_ids_file = channel.value(file("$params.SUBSET_FILE", checkIfExists: true))
        target_bed_files = channel.fromFilePairs("$params.TARGET_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }
        subset_bed_files = subsetBED(target_bed_files, subset_ids_file).subset_bed_files
        
        pcs_and_genotypes = LocoPCA.out.confounders.join(subset_bed_files, failOnDuplicate: true)
        analysis_bed_triplets = subset_bed_files
    } else if (params.DENSE_MAPPING_FILE != "NO_DENSE_MAPPING_FILE") {
        prioritized_variants = channel.value(file("$params.DENSE_MAPPING_FILE", checkIfExists: true))
        imputed_bgen_files = channel.fromFilePairs("$params.BGEN_FILES", size: 3, checkIfExists: true){ f -> f.name.replaceAll(/\.bgen(\.bgi)?$|\.sample$/,'') }
        dense_bed_files = denseBED(imputed_bgen_files, prioritized_variants).dense_bed_files

        pcs_and_genotypes = LocoPCA.out.confounders.join(dense_bed_files, failOnDuplicate: true)
        analysis_bed_triplets = dense_bed_files
    } else {
        pcs_and_genotypes = LocoPCA.out.confounders.join(bed_files, failOnDuplicate: true)
        analysis_bed_triplets = bed_files
    }

    EstimationInputs(
        pcs_and_genotypes,
        LocoPCA.out.traits,
        estimands_file
    )

    // Estimation
    EstimationWorkflow(
        EstimationInputs.out.transpose(), estimator_config
    )

    // TarGWAS Report (HTML + per-estimator summary CSVs)
    if (params.REPORT == true) {
        // Collect all per-chromosome analysis BED triplets, merging when
        // more than one chromosome is present. The subset/dense paths
        // emit a single triplet and bypass the merge.
        analysis_bed_collected = analysis_bed_triplets
            .map { _id, files -> files }
            .collect()
        bed_branched = analysis_bed_collected.branch { all_files ->
            single:   all_files.size() == 3
            multiple: all_files.size() > 3
        }
        merged_multi = mergeBEDS(
            bed_branched.multiple.map { all_files -> tuple('analysis', all_files) }
        ).map { _id, files -> files }
        report_bed = bed_branched.single.mix(merged_multi)

        // Phenotype is optional: only piped in when REPORT_OUTCOME_COL is set.
        report_pheno = params.REPORT_OUTCOME_COL != "NO_REPORT_OUTCOME_COL" ?
            LocoPCA.out.traits :
            channel.value(file("${projectDir}/assets/NO_PHENO"))

        ReportWorkflow(
            EstimationWorkflow.out.merged_hdf5,
            report_bed,
            report_pheno,
        )
    }

    // Generate sieve variance plateau estimates
    if (params.SVP == true) {
        if (params.PREVALENCE != "NO_SET_PREVALENCE") {
            error "SVP is not compatible with a set PREVALENCE parameter."
        }
        // IID Genotypes for SVP (needs all chromosomes merged)
        qc_file = channel.value(file("$params.QC_FILE", checkIfExists: true))
        flashpca_excl_reg = channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS", checkIfExists: true))
        ld_blocks = channel.value(file("$params.LD_BLOCKS", checkIfExists: true))
        IIDGenotypes(
            flashpca_excl_reg,
            ld_blocks,
            bed_files,
            qc_file,
            LocoPCA.out.traits,
        )
        genotypes = IIDGenotypes.out.map{genotypes_id, genotypes -> genotypes}.collect()
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result.collect(), 
            genotypes,
        )
    }
}
