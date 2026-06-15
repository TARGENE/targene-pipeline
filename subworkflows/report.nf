include { TarGenevizReport } from '../modules/report.nf'

workflow ReportWorkflow {
    take:
        hdf5_results        // path  — merged results.hdf5 emitted by GenerateOutputs
        bed_files           // path  — already-merged [bed,bim,fam] triplet or
                            //         single sentinel `NO_BED_PREFIX` file
        pheno_file          // path  — phenotype TSV/CSV or sentinel `NO_PHENO`

    main:
        TarGenevizReport(hdf5_results, bed_files, pheno_file)

    emit:
        report    = TarGenevizReport.out.report
        tabulated = TarGenevizReport.out.tabulated
}
