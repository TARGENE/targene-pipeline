process TarGenevizReport {
    publishDir "${params.OUTDIR}/report", mode: 'copy'
    label 'targeneviz_image'

    input:
        path hdf5
        path bed_files
        path pheno

    output:
        path "report.html",           emit: report
        path "${stem}_tabulated",     emit: tabulated

    script:
        stem       = hdf5.getBaseName()
        forest_top = params.REPORT_FOREST_TOP
        ld_top     = params.REPORT_LD_TOP
        ld_window  = params.REPORT_LD_WINDOW_BP
        bed_list   = bed_files instanceof List ? bed_files : [bed_files]
        has_bed    = !(bed_list.size() == 1 && bed_list[0].getName() == 'NO_BED_PREFIX')
        has_pheno  = pheno.getName() != 'NO_PHENO'
        bed_file   = has_bed ? bed_list.find { it.getName().endsWith('.bed') } : null
        bed_prefix = bed_file ? bed_file.toString().take(bed_file.toString().length() - 4) : ''
        bed_arg    = has_bed   ? "--bed=${bed_prefix}" : ''
        pheno_arg  = has_pheno ? "--pheno=${pheno} --sample-id-col=${params.REPORT_SAMPLE_ID_COL} --outcome-col=${params.REPORT_OUTCOME_COL}" : ''
        """
        TEMPD=\$(mktemp -d)
        # /opt is the precompiled depot; it is made world-writable in the
        # image build so Julia can create cachefile locks at load time.
        # Prepend a per-task scratch depot so any just-in-time compilation
        # output is written there rather than mutating the shared depot.
        JULIA_DEPOT_PATH=\$TEMPD:/opt julia --project=/TarGWAS --startup-file=no \
            /TarGWAS/bin/targeneviz.jl report \
            ${hdf5} \
            report.html \
            ${bed_arg} \
            ${pheno_arg} \
            --ld-window-bp=${ld_window} \
            --ld-top=${ld_top} \
            --forest-top=${forest_top}
        """
}
