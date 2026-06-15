include { PCA } from './pca.nf'
include { EstimationInputs } from '../modules/estimation_inputs.nf'
include { EstimationWorkflow } from '../subworkflows/estimation.nf'
include { SVPWorkflow } from '../subworkflows/svp.nf'
include { ReportWorkflow } from '../subworkflows/report.nf'

workflow TARGENE {
    // Define Parameters
    bgen_files = channel.fromPath("$params.BGEN_FILES", checkIfExists: true).collect().toList()
    estimands_file = channel.value(file("$params.ESTIMANDS_CONFIG"))
    estimator_config = channel.fromPath(EstimatorsConfig.create(params.ESTIMATORS_CONFIG, params.OUTDIR))

    // PCA
    PCA()

    // Estimation Inputs
    pcs_and_genotypes = PCA.out.pcs.combine(bgen_files)
    EstimationInputs(
        pcs_and_genotypes,
        PCA.out.traits,
        estimands_file
    )

    // generate estimates
    EstimationWorkflow(
        EstimationInputs.out.transpose(),
        estimator_config,
    )

    // TarGWAS Report (HTML + per-estimator summary CSVs)
    // TARGENE uses BGEN inputs, so no BED can be piped in. Phenotype is
    // optional: piped in only when REPORT_OUTCOME_COL is set.
    if (params.REPORT == true) {
        report_bed = channel.value(file("${projectDir}/assets/NO_BED_PREFIX"))
        report_pheno = params.REPORT_OUTCOME_COL != "NO_REPORT_OUTCOME_COL" ?
            PCA.out.traits :
            channel.value(file("${projectDir}/assets/NO_PHENO"))
        ReportWorkflow(
            EstimationWorkflow.out.merged_hdf5,
            report_bed,
            report_pheno,
        )
    }

    // Generate sieve variance plateau estimates
    genotypes = PCA.out.iid_genotypes.map{genotypes_id, genotypes -> genotypes}.collect()
    if (params.SVP == true){
        sieve_results = SVPWorkflow(
            EstimationWorkflow.out.hdf5_result.collect(), 
            genotypes,
        )
    }
}
