# GWAS

Perhaps the most common study design in population genetics is the Genome-Wide Association Study (GWAS). In TarGene, a GWAS can only be performed across [plink BED files](https://zzz.bwh.harvard.edu/plink/binary.shtml). A minimalist run configuration for a GWAS in TarGene looks like the following:

```conf
params {
    ESTIMANDS_CONFIG = "gwas_config.yaml"
    ESTIMATORS_CONFIG = "tmle-ose--glm"

    // UK-Biobank specific parameters
    BED_FILES = "unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    UKB_CONFIG = "ukbconfig_gwas.yaml"
    TRAITS_DATASET = "dataset.csv"
    UKB_WITHDRAWAL_LIST = "withdrawal_list.txt"
}
```

All UK-Biobank specific parameters have been described in the [Setup](@ref) section. The only new parameters are the `ESTIMANDS_CONFIG` and the `ESTIMATORS_CONFIG`. These parameters describe the estimands (questions of interest) and how to estimate them respectively. Since we are performing a GWAS, we are interested in the effect of all variants, across all phenotypes in the `UKB_CONFIG` file. Because a GWAS is expensive, this file will only contain 3 traits/phenotypes here.

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

In this case, `Body mass index (BMI)` is the outcome of interest and the `Number of vehicles in household` and `Cheese intake` traits will be used as extra predictors. The reason why this is the case is because they are described as such in the `ESTIMANDS_CONFIG`

```yaml
type: gwas

outcome_extra_covariates:
  - "Number of vehicles in household"
  - "Cheese intake"
```

Where the optional `outcome_extra_covariates` are variables to be used as extra predictors of the outcome (but not as confounders). Importantly, if this section had been omitted, both `Number of vehicles in household` and `Cheese intake` would be considered as outcomes by TarGene, and 3 GWASs instead of 1 would be run.

Finally, the `ESTIMATORS_CONFIG = "tmle-ose--glm"` says that two different estimation strategies will be used. One using the Targeted Minimum-Loss Estimator and one using the One-Step Estimator. The estimtion results can then be compared. Both methods will use a Generalized Linear Model (glm) to learn the outcome models and propensity scores.

The GWAS can be run as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.11.0 -profile local
```

The successful execution will result in the creation of a `results` folder containing at least 3 files:

- results.hdf5
- results.json
- QQ.png

Both `results.hdf5` and `results.json` contain the same estimation results in HDF5 and JSON formats respectively. The `QQ.png` is a Q-Q plot of your results that looks like the following

![GWAS_QQ](../assets/gwas_QQ.png)