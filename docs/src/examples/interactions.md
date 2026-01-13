# Targeted Interaction Study

!!! note "Read First"
    Make sure you have read the [Setup (Read First)](@ref) section.

In TarGene you can investigate the interacting effects between genetic variants and between genetic variants and the environment up to any order (GxG, GxE, GxGxG, ...). However, the computational complexity associated with interaction testing grows exponentially. Interaction studies are thus very often targeted, for instance towards a specific transcription factor `TF`. For this transcription factor, multiple variants may have been identified because they are believed to play different roles. In TarGene such a scenario can be encoded via the `group` mode of the `ESTIMANDS_CONFIG` file. An example is more informative than a long description. Let's save the following content to a `new_interaction_config.yaml` file.

```yaml
type: groups

estimands:
  - type: AIE
    orders: [2, 3]

variants:
  TF:
    variants_1:
      - "1:238411180:T:C"
      - "3:3502414:T:C"
    variants_2:
      - "2:14983:G:A"

extra_treatments:
  - "Cheese intake"

outcome_extra_covariates:
  - "Number of vehicles in household"
```

The `estimands` section describes the types of estimands we are interested in, here interactions are defined by the Average Interaction Effects (`AIE`). We would like to scan both orders 2 and 3.

The `variants` section contains groups of variants for which the interactions will be generated independently. Here there is only one group called `TF`, and it contains two subgroups: `variants_1` and `variants_2`. This means that combinations of `variants_1` and `variants_2` will define the interactions under investigation. That is;  (1:238411180:T:C, 2:14983:G:A) and (3:3502414:T:C, 2:14983:G:A). Why did we set the interaction `orders` to be both 2 and 3 since only pairwise interactions can be generated here?

This is because the `extra_treatments` section can list additional environmental variables, thus increasing the potential interaction order. In this case, there is one: Cheese intake. This means that additional pairwise interactions will be generated: (1:238411180:T:C, Cheese intake), (2:14983:G:A, Cheese intake) and (3:3502414:T:C, Cheese intake). But also the following order 3 points interactions: (1:238411180:T:C, 2:14983:G:A, Cheese intake) and (3:3502414:T:C, 2:14983:G:A, Cheese intake).

!!! note "Note on variant IDs"
    The variant IDs must match the IDs of your BGEN files. Here they are identified by chr:pos:ref:alt as but in your case it may be via the rsID.

Finally, the `outcome_extra_covariates` describes additional variables predictive of the outcome, but that do not confound the effect of treatment variables (`variants` and `extra_treatments`).

Importantly, all traits variables in the `TRAITS_DATASET` that are not listed as either `outcome_extra_covariates` or `extra_treatments` will be considered as outcome variables. The `TRAITS_DATASET` is a plain CSV file that contains an `eid` column that must match the `IID` column in both the `BED_FILES` and `BGEN_FILES` to match individuals.

Finally, let's bring it all together and create the configuration file for Nextflow in our `my_targeted_interaction_study.config` file:

```conf
params {
    ESTIMANDS_CONFIG = "new_interaction_config.yaml"
    BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    BGEN_FILES = "test/assets/unphased_bgen/ukb_chr{1,2,3}.{bgen,bgen.bgi,sample}"
    UKB_CONFIG = "test/assets/ukbconfig_small.yaml"
    TRAITS_DATASET = "test/assets/dataset.csv"
    UKB_WITHDRAWAL_LIST = "test/assets/withdrawal_list.txt"
}
```

As usual, the pipeline can then be run as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.13.0 -profile local -c my_targeted_interaction_study.config -resume
```

When the pipeline terminates, in the `results` directory you should see a `QQ.png` (QQ plot of your results) and a `results.summary.yaml` file of summary statistics. The summary statistics contain entries corresponding to each estimated interaction and of the following form:

```yaml
- EFFECT_TYPE: "AIE"
  OUTCOME: Pork intake
  TREATMENTS:
    - 1:238411180:T:C: "CT => CC"
      Cheese intake: "3 => 2"
    - 1:238411180:T:C: "CT => CC"
      Cheese intake: "2 => 1"
    - 1:238411180:T:C: "CT => CC"
      Cheese intake: "1 => 4"
  WTMLE_TUNEDXGBOOST_TUNEDXGBOOST:
    PVALUE: 0.0009868797600783511
    COMPONENTS:
      - PVALUE: 0.8719094198750913
        EFFECT_SIZE: 0.004502950838730441
      - PVALUE: 0.01141555629312441
        EFFECT_SIZE: 0.06155196300651498
      - PVALUE: 9.74820724988353e-5
        EFFECT_SIZE: -0.10293580652169768
  OSE_TUNEDXGBOOST_TUNEDXGBOOST:
    PVALUE: 0.9268263087987019
    COMPONENTS:
      - PVALUE: 0.9602243491598036
        EFFECT_SIZE: 0.0013928294158895753
      - PVALUE: 0.6078313272814673
        EFFECT_SIZE: 0.01246817359406354
      - PVALUE: 0.5583152167851629
        EFFECT_SIZE: -0.015373244259236892
```

The `TREATMENTS` field contains all variant and environmental variables changes defining the interaction. One such change is the allelic change "CT => CC" for 1:238411180:T:C combined with a Cheese intake change of "2 => 1". Because 1:238411180:T:C and Cheese intake are categorical variables, there are multiple such changes and the full interaction is itself multidimensional. However, if the frequency of a given treatment regimen is too small (`POSITIVITY_CONSTRAINT`, by default 0.01), some components of an interaction may not be estimated. This is the case here since the "TT" genotype is not part of the interaction.

The next sections we are interested in are the actual estimates for those interactions. By default TarGene estimates them using a weighted TMLE estimator and a One Step estimator, both using XGBoost. Those estimators have the same asymptotic behaviour but might differ in finite samples. Since computing the One Step estimate is cheap once the TMLE has been estimated we also report it by default. The corresponding sections are `WTMLE_TUNEDXGBOOST_TUNEDXGBOOST` and `OSE_TUNEDXGBOOST_TUNEDXGBOOST` respectively. First a p-value corresponding to the multidimensional interaction is reported, it corresponds to testing whether the full interaction is non zero. Then each component corresponding to the different treatment changes are reported in the same order with both a p-value and an effect size.  We can see that the two estimators indeed provide very different p-values here in this small and artificial dataset. In the UK-Biobank, it was shown that they behave very similarly in most [usual conditions](https://academic.oup.com/biostatistics/article/26/1/kxaf030/8268272).

For more on TarGene's outputs, see [Understanding TarGene's Outputs](@ref).