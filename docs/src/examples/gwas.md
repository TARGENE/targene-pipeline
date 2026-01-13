# GWAS

!!! note "Read First"
    Make sure you have read the [Setup (Read First)](@ref) section.

Perhaps the most common study design in population genetics is the Genome-Wide Association Study (GWAS). In TarGene, for computational reasons, a GWAS can only be performed across [plink BED files](https://zzz.bwh.harvard.edu/plink/binary.shtml). For our GWAS example, let's save the following snippet to a `my_gwas.conf` file:

```conf
params {
    BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    ESTIMATORS_CONFIG = "wtmle--glm"
    TRAITS_DATASET = "test/assets/dataset.csv"
    UKB_WITHDRAWAL_LIST = "test/assets/withdrawal_list.txt"
    ESTIMANDS_CONFIG = "test/assets/gwas_config.yaml"
    UKB_CONFIG = "test/assets/ukbconfig_gwas.yaml"
}
```

Let's now unpack what this configuration file is saying. Since we have already described most parameters in the [Setup (Read First)](@ref) section, we focus on the two remaining parameters:

- `ESTIMATORS_CONFIG`: Defines the estimation strategy. We will be using Targeted Minimum-Loss Estimator with a simple Generalized Linear Model (glm) to learn the outcome models (``Q_Y``)and propensity scores (``G``). For small studies, more complex machine-learning models can be used but for a large GWAS this might be computationally too intensive. See [Specifying a Targeted Estimator](@ref) for more information on this.
- `ESTIMANDS_CONFIG`: Since we are performing a GWAS, we are interested in the effect of all variants, across all phenotypes in the `UKB_CONFIG` file.  The `ESTIMANDS_CONFIG` is a pretty succinct YAML file for a GWAS, it could contain only one line. Here we will be a little more fancy and add extra outcome predictors to improve the precision of inference. The `test/assets/gwas_config.yaml` contains the following which is self explanatory.

```yaml
type: gwas

outcome_extra_covariates:
  - "Number of vehicles in household"
  - "Cheese intake"
```

The optional `outcome_extra_covariates` are variables to be used as extra predictors of the outcome (but not as confounders). In traditional GWAS these are typically age and sex.

The GWAS can then be run as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.13.0 -profile local -c my_gwas.config -resume
```

The workflow should start running and download the first container image. Once successfully completed, the workflow results will be found in the `results` folder. The main output is the summary stats YAML file (`results/results.summary.yaml`) which contains entries of the form:

```yaml
- EFFECT_TYPE: "ATE"
  OUTCOME: Body mass index (BMI)
  TREATMENTS:
    - 1:235339932:G:A: "0 => 1"
    - 1:235339932:G:A: "1 => 2"
  WTMLE_GLM_GLM:
    PVALUE: 0.67290478200139
    COMPONENTS:
      - PVALUE: 0.44667770202884105
        EFFECT_SIZE: -0.5283205242017841
      - PVALUE: 0.5038425323927094
        EFFECT_SIZE: 0.3319580826783738
```

Remember that in TarGene we estimate the effect of allelic changes and we do not assume that the effect of adding one minor allele is the same as the effect of adding the second minor allele. Hence each variant as up to two distinct effect sizes each with their own p-value. The overall p-value (under the `WTMLE_GLM_GLM` section) comes from a multidimensional test asking whether **any** of the allelic changes is non zero. In some cases, only one of the two changes is estimated. This is because there wasn't sufficient data to be confident about assessing the second allelic change (see the `POSITIVITY_CONSTRAINT` in the workflow parameters).

For more on TarGene's outputs, see [Understanding TarGene's Outputs](@ref).