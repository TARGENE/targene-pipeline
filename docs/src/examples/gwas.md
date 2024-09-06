# GWAS

Perhaps the most common study design in population genetics is the Genome-Wide Association Study (GWAS). In TarGene, a GWAS can only be performed across [plink BED files](https://zzz.bwh.harvard.edu/plink/binary.shtml). A minimalist run configuration for a GWAS in TarGene looks like the following:

```conf
params {
    ESTIMANDS_CONFIG = "gwas_config.yaml"
    ESTIMATORS_CONFIG = "wtmle--glm"

    // UK-Biobank specific parameters
    BED_FILES = "unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    UKB_CONFIG = "ukbconfig_gwas.yaml"
    TRAITS_DATASET = "dataset.csv"
    UKB_WITHDRAWAL_LIST = "withdrawal_list.txt"
}
```

Apart from the data related parameters, there are two main parameters here: `ESTIMANDS_CONFIG` and the `ESTIMATORS_CONFIG`. These parameters describe the estimands (questions of interest) and how to estimate them respectively. Since we are performing a GWAS, we are interested in the effect of all variants, across all phenotypes in the `UKB_CONFIG` file. 

The `ESTIMANDS_CONFIG` is a pretty succinct YAML file for a GWAS, it could contain only one line. Here we will be a little more fancy and add extra outcome predictors to improve the precision of inference

```yaml
type: gwas

outcome_extra_covariates:
  - "Number of vehicles in household"
  - "Cheese intake"
```

The optional `outcome_extra_covariates` are variables to be used as extra predictors of the outcome (but not as confounders). How does TarGene know how to extract these traits from the the UK Biobank main dataset? This is thanks to the `UKB_CONFIG` file, which maps UK Biobank data fields to traits (see [The `UKB_CONFIG` Configuration File](@ref)). In our case, this file will contain both the outcome of interest (BMI) and the extra predictors we need (Number of vehicles in household and Cheese intake).

```yaml
traits:
  - fields:
      - "21001"
    phenotypes:
      - name: "Body mass index (BMI)"
  - fields:
      - "728"
    phenotypes:
      - name: "Number of vehicles in household"
  - fields:
      - "1408"
    phenotypes:
      - name: "Cheese intake"
```

Note that the outcome is defined implicitely, any trait in the `UKB_CONFIG` file which is not in the `outcome_extra_covariates` will be considered as an outcome. So you can run multiple GWAS at once by simply adding another trait definition to the above file.

Finally, the `ESTIMATORS_CONFIG = "wtmle--glm"` defines the estimation strategy. We will be using Targeted Minimum-Loss Estimator with a simple Generalized Linear Model (glm) to learn the outcome models (``Q_Y``)and propensity scores (``G``).

The GWAS can then be run as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.11.0 -profile local
```