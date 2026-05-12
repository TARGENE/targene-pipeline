def filepath_matches_chr_prefix(fp, chr_prefix){
    def fp_string = fp.normalize().toString()
    return fp_string.contains(chr_prefix.normalize().toString() + ".")
}

def leave_chr_out(chr_prefix, bed_files){
    def bed_files_not_matching_chr_prefix = bed_files.findAll{ fp -> !filepath_matches_chr_prefix(fp, chr_prefix) }
    return [chr_prefix, bed_files_not_matching_chr_prefix]
}
