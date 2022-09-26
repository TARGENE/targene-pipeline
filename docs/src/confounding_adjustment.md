# Confounding Adjustment

To account for potential confounding effect due to population stratification, we extract principal components from the genetic data using [flashpca](https://github.com/gabraham/flashpca). We follow the recommended procedure for this tool which implies some preprocessing and filtering. The following arguments are compulsory:

- `LD_BLOCKS`: A path to pre-identified linkage desequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data.
- `FLASHPCA_EXCLUSION_REGIONS`: A path to the flashpca special exclusion regions which is provided in their repository.
- `NB_PCS` (default: 6): The number of PCA components to extract.
