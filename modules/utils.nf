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

def processEstimatorsConfig(configValue) {
    // If it's a list, take the first element (ex: in NULL simulation test config)
    if (configValue instanceof List) {
        if (configValue.size() == 0) {
            error "ESTIMATORS_CONFIG list is empty"
        }
        configValue = configValue[0]
    }

    // At this point, configValue should be a string
    if (!(configValue instanceof String)) {
        error "ESTIMATORS_CONFIG must be a string, a list containing a string, or a file path"
    }

    // Remove brackets and quotes if present
    configValue = configValue.replaceAll(/^\[|\]$|"/, '')

    def configFile = file(configValue)

    // If it's not an existing file, create an empty file with this name
    if (!configFile.exists()) {
        configFile = file("${launchDir}/${configValue}")
        configFile.text = '' // Create an empty file
        println "Created empty estimators config file: ${configFile}"
    } else {
        println "Using existing estimators config file: ${configFile}"
    }

    return configFile.toString()
}