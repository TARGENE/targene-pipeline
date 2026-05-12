class LongestPrefix {
    static compute(files){
        // Only one file, strangely it is not passed as a list
        if (files instanceof Collection == false) {
            return files.getName()
        }
        // More than one file
        def index = 0
        while(true){
            def current_prefix = files[0].getName()[0..index]
            for (file in files){
                if(file.getName()[0..index] != current_prefix){
                    return current_prefix[0..-2]
                }
            }
            index++
        }
    }
}

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