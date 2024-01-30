include { longest_prefix } from './utils.nf'

process MakeDataset {
    container "olivierlabayle/tl-core:cvtmle"
    publishDir "${params.OUTDIR}", mode: 'symlink'
    label "bigmem"

    input:
        path bgenfiles
        path traits
        path confounders
    
    output:
        path "dataset.arrow"
    
    script:
        prefix = longest_prefix(bgenfiles)
        """
        """
}