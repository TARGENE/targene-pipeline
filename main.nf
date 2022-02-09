#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.PHENOTYPES_LIST = "NONE"
params.QUERIES_MODE = "given"
params.PHENOTYPES_IN_PARALLEL = false
params.THRESHOLD = 0.9
params.CROSSVAL = false

include { IIDGenotypes } from './modules/genotypes.nf'
include { generatePCs } from './modules/confounders.nf'
include { queriesFromASBxTransActors; filterASB; queriesFromQueryFiles } from './modules/queries.nf'
include { phenotypesFromGeneAtlas } from './modules/phenotypes.nf'
include { VariantRun as TMLE; VariantRun as CrossVal} from './modules/tmle.nf'
include { GRMPart; AggregateIDFiles; PrependSize } from './modules/grm.nf'


def grm_criteria = branchCriteria {
                id: it.getName().endsWith(".grm.id")
                bin: it.getName().endsWith(".grm.bin")
                n_bin: it.getName().endsWith(".grm.N.bin")
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

workflow generateGRM {
    take:
        iid_genotypes

    main:
        grm_parts = Channel.from( 1..params.GRM_NSPLITS )
        GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)

        GRMPart.flatten().branch {
            id: it.getName().endsWith(".grm.id")
            bin: it.getName().endsWith(".grm.bin")
            n_bin: it.getName().endsWith(".grm.N.bin")
        }
        .set { result }

        // Aggregate ID files
        AggregateIDFiles(result.id.collect())

        // Prepend size to bin files
        PrependSize(result.bin)

    emit:
        grm_ids = AggregateIDFiles.out
        grm_bin = PrependSize.out

}


workflow generateConfounders {
    take:
        iid_genotypes

    main:
        generatePCs(iid_genotypes.collect())

    emit:
        generatePCs.out
}


workflow generateQueries{
    if (params.QUERIES_MODE == "ASBxTransActors") {
        asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
        trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
        excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
        bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
        filterASB(asb_snp_ch.collect())
        queries = queriesFromASBxTransActors(filterASB.out.filtered_asb_snps, 
                                             trans_actors, 
                                             bgen_files_ch.toList(), 
                                             excluded_snps)
    }
    else if (params.QUERIES_MODE == "given"){
        queries = queriesFromQueryFiles(Channel.fromPath("$params.QUERY_FILES", checkIfExists: true))
    }
    emit:
        queries

}


workflow generatePhenotypes {
    binary_phenotypes = Channel.value(file("$params.BINARY_PHENOTYPES"))
    continuous_phenotypes = Channel.value(file("$params.CONTINUOUS_PHENOTYPES"))
    bridge = Channel.value(file("$params.GENEATLAS_BRIDGE"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))

    phenotypesFromGeneAtlas(binary_phenotypes, continuous_phenotypes, bridge, withdrawal_list)
    
    emit:
        phenotypesFromGeneAtlas.out
}


workflow generateEstimates {
    take:
        phenotypes_file
        queries_files
        confounders_file

    main:
        bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))

        if (params.PHENOTYPES_LIST != "NONE" & params.PHENOTYPES_IN_PARALLEL == true){
            phenotypes_list = Channel.fromPath("$params.PHENOTYPES_LIST", checkIfExists: true)
                                    .splitText(file: true)
        }
        else if (params.PHENOTYPES_LIST != "NONE" & params.PHENOTYPES_IN_PARALLEL == false) {
            phenotypes_list = Channel.fromPath("$params.PHENOTYPES_LIST", checkIfExists: true)
        }
        else if (params.PHENOTYPES_LIST == "NONE") {
            phenotypes_list = phenotypes_file.splitCsv(sep: ",", limit: 1)
                                    .flatten()
                                    .filter { it != "eid"}

            if (params.PHENOTYPES_IN_PARALLEL == false) {
                phenotypes_list = phenotypes_list.collectFile(name: 'phenotypes_list.txt', newLine: true)
            }
            else {
                count = 0
                phenotypes_list = phenotypes_list
                        .collectFile(){ item -> [ "phenotypes_${count++}.txt", item]}
            }
        }

        phen_list_to_queries = queries_files.combine(phenotypes_list)
        // compute TMLE estimates
        TMLE("epistasis", bgen_files_ch.collect(), phenotypes_file, confounders_file, estimator_file, phen_list_to_queries)
        // Aggregate results
        TMLE.out.collectFile(name:"epistasis_estimates.csv",
                            keepHeader: true,
                            skip: 1,
                            storeDir: "$params.OUTDIR")
        
        if (params.CROSSVAL == true) {
            CrossVal("crossval", bgen_files_ch.collect(), phenotypes_file, confounders_file, estimator_file, phen_list_to_queries)
            // Aggregate results
            CrossVal.out.collectFile(name:"crossval_estimates.csv",
                                keepHeader: true,
                                skip: 1,
                                storeDir: "$params.OUTDIR")
        }


}


workflow {
    // Generate queries
    generateQueries()

    // Generate IID Genotypes
    generateIIDGenotypes()

    // generate confounders
    // generateConfounders(generateIIDGenotypes.out)

    // generate GRM
    generateGRM(generateIIDGenotypes.out)

    // generate phenotypes
    // generatePhenotypes()

    // // generate estimates
    // generateEstimates(generatePhenotypes.out, generateQueries.out.flatten(), generateConfounders.out)
    
}