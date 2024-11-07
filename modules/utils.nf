def longest_prefix(files){
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

def filepath_matches_chr_prefix(fp, chr_prefix){
    def fp_string = fp.normalize().toString()
    return fp_string.contains(chr_prefix.normalize().toString() + ".")
}

def leave_chr_out(chr_prefix, bed_files){
    def bed_files_not_matching_chr_prefix = bed_files.findAll{ fp -> !filepath_matches_chr_prefix(fp, chr_prefix) }
    return [chr_prefix, bed_files_not_matching_chr_prefix]
}

def CreateEstimatorsConfigChannel(configValue) {
    configFiles = []
    // Ensure configValue is a list
    if (!(configValue instanceof List)) {
        configValue = [configValue]
    } 

    // Return the list of estimator names or paths
    for (estimator in configValue) {
        def configFile = file(estimator)
        // If it's not an existing file, create an empty file with this name
        if (!configFile.exists()) {
            configFile = file("${params.OUTDIR}/${estimator}")
            configFile.text = '' // Create an empty file
        }
        configFiles.add(configFile)
    }

    return Channel.fromList(configFiles).collect()
/*
    return Channel.fromList(configValueList)
        .map { estimator ->
            def configFile = file(estimator)
            if (!configFile.exists()) {
                configFile = file("${params.OUTDIR}/${estimator}")
                configFile.text = '' // Create an empty file
            }
            return configFile
        }
        .collect() */
}

/*
def processEstimatorsConfig(configValue) {
    def configFiles = []

    // Ensure configValue is a list
    if (!(configValue instanceof List)) {
        configValue = [configValue]
    }

    for (estimator in configValue) {
        def configFile = file(estimator)

        // If it's not an existing file, create an empty file with this name
        if (!configFile.exists()) {
            configFile = file("${params.OUTDIR}/${estimator}")
            configFile.text = '' // Create an empty file
        }

        configFiles.add(configFile)
    }

    // Return the list of configFiles
    println configFiles
    return configFiles
}
*/
/*
def processEstimatorsConfig(configValue) {
    if (configValue instanceof List) {
        if (configValue.size() == 0) {
            error "ESTIMATORS_CONFIG list is empty"
        }
        configValue = configValue[0]
    }

    def configFile = file(configValue)

    // If it's not an existing file, create an empty file with this name
    if (!configFile.exists()) {
        configFile = file("${launchDir}/${configValue}")
        configFile.text = '' // Create an empty file
    } 

    return configFile.toString()
}
*/