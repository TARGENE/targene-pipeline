# PCA Adjustment

To account for potential confounding effect due to population stratification, we extract principal components from the genetic data using [flashpca](https://github.com/gabraham/flashpca). We follow the recommended procedure for this tool which implies some preprocessing and filtering.

In principle, you shouldn't need to change these parameters apart from the number of principal components (`NB_PCS`).

- `NB_PCS` (optional, default: 6): The number of PCA components to extract.
- `LD_BLOCKS` (optional): A path to pre-identified linkage disequilibrium blocks around the variants that will be queried for causal effect estimation. Those LD blocks will be removed from the data used for PCA.
- `FLASHPCA_EXCLUSION_REGIONS` (required, default: [assets/exclusion_regions_hg19.txt](https://github.com/TARGENE/targene-pipeline/assets/exclusion_regions_hg19.txt)): A path to the flashpca special exclusion regions which is provided in [their repository](https://github.com/gabraham/flashpca/blob/master/exclusion_regions_hg19.txt).
- `MAF_THRESHOLD` (optional, default: 0.01): Only variants with that minor allele frequency are used to compute principal components.
