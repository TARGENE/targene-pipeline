#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.FIELDS_METADATA = "NO_FILE"
params.EXTRA_COVARIATES = "NO_FILE"
params.ENCRYPTION_FILE = "NO_FILE"
params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.PHENOTYPES_LIST = "NO_FILE"
params.QUERIES_MODE = "given"
params.CALL_THRESHOLD = 0.9
params.MINOR_CAT_FREQUENCY = 0.001
params.SAVE_FULL = false
params.PHENOTYPES_BATCH_SIZE = 1
params.GRM_NSPLITS = 100
params.MAF_THRESHOLD = 0.01
params.NB_PCS = 6
params.NB_VAR_ESTIMATORS = 0
params.MAX_TAU = 0.8
params.PVAL_SIEVE = 0.05
params.OUTDIR = "$launchDir/results"

include { IIDGenotypes } from './modules/genotypes.nf'
include { generatePCs; MergeExtraCovariatesAndPCs } from './modules/covariates.nf'
include { FromASBxTransActors; FromGivenQueries } from './modules/queries.nf'
include { phenotypesFromGeneAtlas as BridgeContinuous; phenotypesFromGeneAtlas as BridgeBinary } from './modules/phenotypes.nf'
include { TMLE as TMLEContinuous; TMLE as TMLEBinary} from './modules/tmle.nf'
include { PhenotypesBatches as ContinuousPhenotypesBatches; PhenotypesBatches as BinaryPhenotypesBatches} from './modules/tmle.nf'
include { GRMPart; AggregateGRM } from './modules/grm.nf'
include { SieveVarianceEstimation } from './modules/sieve_variance.nf'
include { Summary } from './modules/summary.nf'


def NbPhenotypes() {
    if (params.PHENOTYPES_LIST != "NO_FILE") {
        reader = file(params.PHENOTYPES_LIST).newReader()
        int lines = 0
        while (reader.readLine() != null) { 
            lines++
        }
        return lines
    }
    else {
        binReader = file(params.BINARY_PHENOTYPES).newReader()
        nbBin = binReader.readLine().split(" ").size()

        contReader = file(params.CONTINUOUS_PHENOTYPES).newReader()
        nbCont = contReader.readLine().split(" ").size()
        // Remove twice the FID and IID columns
        return nbBin + nbCont - 4
    }
}

NB_PHENOTYPES = NbPhenotypes()

process UKBConv {
    container "olivierlabayle/ukbmain:v0.2.0"

    input:
        path fields_file
        path encr_file

    output:
        path "output.csv"

    script:
        """
        ukbconv $encr_file csv -i$fields_file -ooutput
        """
}

process DecodeMainDataset {
    container "olivierlabayle/ukbmain:v0.2.0"

    input:
        path dataset_file
        path fields_metadata

    output:
        path "main_dataset.csv"

    script:
        """
        julia --project=/UKBMain/ --startup-file=no /UKBMain/bin/decode.jl \
        $dataset_file $fields_metadata main_dataset.csv
        """
}


workflow generateIIDGenotypes {
    qc_file = Channel.value(file("$params.QC_FILE"))
    flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
    ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
    bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

    IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file)

    emit:
        IIDGenotypes.out
}


workflow generateCovariates {
    take:
        iid_genotypes

    main:
        covariates = generatePCs(iid_genotypes.collect())
        if (params.EXTRA_COVARIATES != "NO_FILE") {
            covariates_file = Channel.value(file("$params.EXTRA_COVARIATES"))
            encryption_file = Channel.value(file("$params.ENCRYPTION_FILE"))
            fields_metadata = Channel.value(file("$params.FIELDS_METADATA"))
            UKBConv(covariates_file, encryption_file)
            DecodeMainDataset(UKBConv.out, fields_metadata)
            covariates = MergeExtraCovariatesAndPCs(covariates, DecodeMainDataset.out)
        }

    emit:
        covariates
}


workflow generateQueriesAndGenotypes{
    bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
    excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
    if (params.QUERIES_MODE == "ASBxTransActors") {
        asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
        outputs = FromASBxTransActors(bgen_files_ch.collect(),
                                             asb_snp_ch.collect(), 
                                             trans_actors, 
                                             excluded_snps)
    }
    else if (params.QUERIES_MODE == "given"){
        query_files = Channel.fromPath("$params.QUERY_FILES", checkIfExists: true).collect()
        outputs = FromGivenQueries(bgen_files_ch.collect(), query_files, excluded_snps)
    }

    emit:
        genotypes = outputs.genotypes
        queries = outputs.queries

}


workflow generatePhenotypes {
    binary_phenotypes = Channel.fromPath("$params.BINARY_PHENOTYPES")
    continuous_phenotypes = Channel.fromPath("$params.CONTINUOUS_PHENOTYPES")
    bridge = Channel.fromPath("$params.GENEATLAS_BRIDGE")
    withdrawal_list = Channel.fromPath("$params.WITHDRAWAL_LIST")
    phenotypes_list = Channel.fromPath("$params.PHENOTYPES_LIST")

    BridgeContinuous(continuous_phenotypes, bridge, withdrawal_list, phenotypes_list, "continuous_phenotypes.csv")

    BridgeBinary(binary_phenotypes, bridge, withdrawal_list, phenotypes_list, "binary_phenotypes.csv")
    
    emit:
        continuous = BridgeContinuous.out
        binary = BridgeBinary.out
}


workflow generateEstimates {
    take:
        genotypes_file
        queries_files
        continuous_phenotypes_file
        binary_phenotypes_file
        covariates_file

    main:
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))

        // compute TMLE estimates for continuous targets
        ContinuousPhenotypesBatches(continuous_phenotypes_file)
        queries_to_continuous_phenotype_batches = queries_files.combine(ContinuousPhenotypesBatches.out.flatten())
        TMLEContinuous(genotypes_file, continuous_phenotypes_file, covariates_file, estimator_file, queries_to_continuous_phenotype_batches, "Real")
        
        // compute TMLE estimates for binary targets
        BinaryPhenotypesBatches(binary_phenotypes_file)
        queries_to_binary_phenotype_batches = queries_files.combine(BinaryPhenotypesBatches.out.flatten())
        TMLEBinary(genotypes_file, binary_phenotypes_file, covariates_file, estimator_file, queries_to_binary_phenotype_batches, "Bool")

        hdf5_files = TMLEContinuous.out.flatten()
                        .concat(TMLEBinary.out.flatten())
                        .map { it -> [it.getName().split("_batch")[0], it]}
                        .groupTuple(size: NB_PHENOTYPES)

    emit:
        hdf5_files
}


workflow generateSieveEstimates {
    take:
        snps_tmle_files
        iid_genotypes
    
    main:
        Channel.from([4, 5, 6]).view()
        if (params.NB_VAR_ESTIMATORS != 0){
            // Build the GRM
            grm_parts = Channel.from( 1..params.GRM_NSPLITS )
            GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)
            AggregateGRM(GRMPart.out.collect())
            // Debug
            Channel.from([7, 8, 9]).view()
            snps_tmle_files.view()
            // Sieve estimation
            sieve_estimates = SieveVarianceEstimation(snps_tmle_files, AggregateGRM.out.grm_ids, AggregateGRM.out.grm_matrix)
        }
        else {
            sieve_estimates = snps_tmle_files.map(it -> [it[0], "NO_FILE"])
        }
    emit:
        sieve_estimates
}

workflow generateSummaries {
    take:
        tmle_files
        sieve_files
    
    main:
        // joining on the prefix which corresponds to a tuple of SNPS
        Summary(tmle_files.join(sieve_files))
}

workflow {
    // Generate queries
    generateQueriesAndGenotypes()

    // Generate IID Genotypes
    generateIIDGenotypes()

    // generate covariates
    generateCovariates(generateIIDGenotypes.out)

    // generate phenotypes
    generatePhenotypes()

    // generate estimates
    generateEstimates(
        generateQueriesAndGenotypes.out.genotypes.first(),
        generateQueriesAndGenotypes.out.queries.flatten(),
        generatePhenotypes.out.continuous.first(),
        generatePhenotypes.out.binary.first(),
        generateCovariates.out.first()
    )

    // generate sieve estimates
    Channel.from([1,2,3]).view()
    generateEstimates.out.view()
    generateSieveEstimates(generateEstimates.out, generateIIDGenotypes.out)

    // generate Summaries
    generateSummaries(generateEstimates.out, generateSieveEstimates.out)
}