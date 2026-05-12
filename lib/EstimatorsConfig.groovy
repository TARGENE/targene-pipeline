class EstimatorsConfig {

    static List create(configValue, outdir) {

        def values = (configValue instanceof List) ? configValue : [configValue]

        def result = []

        for (estimator in values) {

            def configFile = new File(estimator)
            def alreadyCreatedFile = new File("${outdir}/${estimator}")

            if (!configFile.exists() && !alreadyCreatedFile.exists()) {
                new File(outdir).mkdirs()
                def f = new File("${outdir}/${estimator}")
                f.text = ''
                result << f
            } else if (alreadyCreatedFile.exists()) {
                result << alreadyCreatedFile
            } else {
                result << configFile
            }
        }

        return result
    }
}