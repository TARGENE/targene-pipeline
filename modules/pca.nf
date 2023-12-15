workflow PCA {
    // Extract traits
    extractTraits()

    // Generate IID Genotypes
    IIDGenotypes(extractTraits.out)

    // Genetic confounders up to NB_PCS
    geneticConfounders(IIDGenotypes.out) 

}