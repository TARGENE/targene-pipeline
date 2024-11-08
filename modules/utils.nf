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
    estimators_ch = Channel.empty()
    // Ensure configValue is a list
    if (!(configValue instanceof List)) {
        configValue = [configValue]
    } 

    // Iterate through each value of this list
    for (estimator in configValue) {
        def configFile = file(estimator)
        // If it's not an existing file, create an empty file with this name
        if (!configFile.exists()) {
            file(params.OUTDIR).mkdirs() // Create OUTDIR if it doesn't exist
            configFile = file("${params.OUTDIR}/${estimator}")
            configFile.text = '' // Create an empty file
        }
        estimators_ch = estimators_ch.mix(Channel.value(configFile))
    }

    return estimators_ch
}
