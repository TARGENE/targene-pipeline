# GWAS

Perhaps the most common study design in population genetics is the Genome-Wide Association Study (GWAS). In TarGene, a GWAS can only be performed across [plink BED files](https://zzz.bwh.harvard.edu/plink/binary.shtml). For our GWAS example, let's save the following snippet to a `my_gwas.conf` file:

```conf
params {
    ESTIMANDS_CONFIG = "test/assets/gwas_config.yaml"
    ESTIMATORS_CONFIG = "wtmle--glm"

    // UK-Biobank specific parameters
    BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    UKB_CONFIG = "test/assets/ukbconfig_gwas.yaml"
    TRAITS_DATASET = "test/assets/dataset.csv"
    UKB_WITHDRAWAL_LIST = "test/assets/withdrawal_list.txt"
}
```

Let's now unpack what this configuration file is saying.

- `BED_FILES`: This parameter points to the plink bed files containing genetic variants information. Notice the "{1,2,3}.{bed,bim,fam}" suffix of the parameter which is telling Nextflow to group chromosome files together.
- `ESTIMATORS_CONFIG`: Defines the estimation strategy. We will be using Targeted Minimum-Loss Estimator with a simple Generalized Linear Model (glm) to learn the outcome models (``Q_Y``)and propensity scores (``G``). See [Specifying a Targeted Estimator](@ref) for more information on this.
- `TRAITS_DATASET`: This is a CSV file containing both covariates and phenotypes of interest.It should contain an `eid` column to uniquely identify individuals and matching the plink individuals' `IID`.
- `UKB_WITHDRAWAL_LIST`: List of individuals to be removed from the analysis.
- `ESTIMANDS_CONFIG`: Since we are performing a GWAS, we are interested in the effect of all variants, across all phenotypes in the `UKB_CONFIG` file.  The `ESTIMANDS_CONFIG` is a pretty succinct YAML file for a GWAS, it could contain only one line. Here we will be a little more fancy and add extra outcome predictors to improve the precision of inference. The `test/assets/gwas_config.yaml` contains the following which is self explanatory.

```yaml
type: gwas

outcome_extra_covariates:
  - "Number of vehicles in household"
  - "Cheese intake"
```

The optional `outcome_extra_covariates` are variables to be used as extra predictors of the outcome (but not as confounders). 

- `UKB_CONFIG`: How does TarGene know how to extract covariates and phenotypes from the the UK Biobank main dataset? This is thanks to the `UKB_CONFIG` file, which maps UK Biobank data fields to traits (see [The `UKB_CONFIG` Configuration File](@ref)). In our case, this file (`test/assets/ukbconfig_gwas.yaml`) contains both the outcome of interest (BMI) and the extra predictors we need (Number of vehicles in household and Cheese intake).

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

The GWAS can then be run as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.11.1 -profile local -c my_gwas.config
```