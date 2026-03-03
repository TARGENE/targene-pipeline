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
    - 1:235339932:G:A: "GG => GA"
    - 1:235339932:G:A: "GA => AA"
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

## Additional Modes

After running the initial GWAS, certain genomic regions may become of interest for further investigation. However, the sparse coverage of directly genotyped data from BED files does not enable post-GWAS analyses such as fine-mapping and colocalisation, that require dense coverage of the loci of interest. To obtain this dense coverage, one can specify a path to imputed BGEN files in the configuration file. However, due to computational regions, performing estimation across an entire set of imputed variants is not feasible for TarGene. To circumvent this, one can specify the `DENSE_MAPPING_FILE` which is a text file with the following format:

```txt
chromosome position label
1 156365093 CHR1
2 14983 CHR2
3 197645995 CHR3
```
and the `WINDOW_SIZE` parameter which sets for the base range around the loci in the `DENSE_MAPPING_FILE`. An example configuration file would appear as such:
```conf
params {
    BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    BGEN_FILES = "test/assets/unphased_bgen/ukb_chr{1,2,3}.{bgen, bgen.bgi, sample}"
    DENSE_MAPPING_FILE = "test/assets/prioritised_snps.txt"
    WINDOW_SIZE = 500000
    ESTIMATORS_CONFIG = "wtmle--glm"
    TRAITS_DATASET = "test/assets/dataset.csv"
    UKB_WITHDRAWAL_LIST = "test/assets/withdrawal_list.txt"
    ESTIMANDS_CONFIG = "test/assets/gwas_config.yaml"
    UKB_CONFIG = "test/assets/ukbconfig_gwas.yaml"
}
```

Additionally, there is a option to investigate loci of interest under GWAS settings without having to run the entire GWAS again. This can be done by specifying `TARGET_BED_FILES` and `SUBSET_FILE` in the configuration file. The `TARGET_BED_FILES` is the set of BED files in which your variants of interest reside whereas the original `BED_FILES` parameter is the full set now used for the computation of PCs to capture population stratification. The `SUBSET_FILE` is simply a list of variants as they appear in the plink BED format:

```text
1:1593249:G:A
1:1816878:C:A
1:3579379:C:G
2:14983:G:A
2:1786744:T:C
2:5274665:A:G
```

The configuration file appears as such:

```conf
params {
    BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    TARGET_BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2}.{bed,bim,fam}"
    SUBSET_FILE = "test/assets/prioritised_snps.txt"
    ESTIMATORS_CONFIG = "wtmle--glm"
    TRAITS_DATASET = "test/assets/dataset.csv"
    UKB_WITHDRAWAL_LIST = "test/assets/withdrawal_list.txt"
    ESTIMANDS_CONFIG = "test/assets/gwas_config.yaml"
    UKB_CONFIG = "test/assets/ukbconfig_gwas.yaml"
}
```

# GWIS

Provided there is sufficient power, one may concern themselves with interaction effects at genome-wide scale or a genome-wide interaction study (GWIS). Whether this query concerns how variants interact with other variants (epistasis) or gene-by-environment interaction, both can be addressed using the GWAS configuration. All that is required for this that differs from the traditional GWAS configuation is some minor additions to the `ESTIMANDS_CONFIG` yaml file.

We can begin by expanding upon our intial estimands configuration. In the first example, we did not need to specify the `estimands` field as they are by default set to be the single order Average Treatment Effects (ATE). Now, that we want more complex estimands we must define the field `estimands` with the subfields `type` and `orders`. Here, we specify that we want to estimate the Average Interaction Effect (AIE) of the second order between each variant and some treatment variable on outcome. These treatment variables are defined under the `extra_treatments` field where we have specified both the genetic variant `1:235342235:C:T` and the environmental variable `Hot drink temperature`. This configuration will then define the AIE between each genetic variant the levels of `1:235342235:C:T` on outcome as well as the AIE of each genetic variant and the levels of `Hot drink temperature`.

```yaml
type: gwas

estimands:
  - type: AIE
    orders: [2]

extra_treatments:
  - Hot drink temperature
  - 1:235342235:C:T

outcome_extra_covariates:
  - "Number of vehicles in household"
  - "Cheese intake"
```

Futhermore, if we were to add and additional order, say `3`:

```yaml
type: gwas

estimands:
  - type: AIE
    orders: [2,3]

extra_treatments:
  - Hot drink temperature
  - 1:235342235:C:T

outcome_extra_covariates:
  - "Number of vehicles in household"
  - "Cheese intake"
```

We would then have the previously defined estimands as well as a three-point interaction between each variant and the levels of both of the values in `extra_treatments`.

The GWAS can then be run as before:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.13.0 -profile local -c my_gwas.config -resume
```

And the output would be similar, however the `EFFECT_TYPE` and `TREATMENTS` will, of course, be different. Here is an example of output for the two point interaction between two variants:

```yaml
- EFFECT_TYPE: "AIE"
  OUTCOME: Body mass index (BMI)
  TREATMENTS:
    - 1:235339932:G:A: "GG => GA"
      1:235342235:C:T: "CC => CT"
    - 1:235339932:G:A: "GA => AA"
      1:235342235:C:T: "CC => CT"
    - 1:235339932:G:A: "GG => GA"
      1:235342235:C:T: "CT => TT"
    - 1:235339932:G:A: "GA => AA"
      1:235342235:C:T: "CT => TT"
  WTMLE_GLM_GLM:
    PVALUE: 0.19908549213021756
    COMPONENTS:
      - PVALUE: 0.09118230780188791
        EFFECT_SIZE: -0.0032775901441661776
      - PVALUE: 0.26941460025816427
        EFFECT_SIZE: 0.0015468229348041384
      - PVALUE: 0.54123876102934512
        EFFECT_SIZE: 0.0008921456782310245
      - PVALUE: 0.71034521098765432
        EFFECT_SIZE: -0.0004356789012345678
```