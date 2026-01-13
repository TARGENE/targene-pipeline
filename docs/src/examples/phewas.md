# PheWAS

!!! note "Read First"
    Make sure you have read the [Setup (Read First)](@ref) section.

A phenome-wide association study (PheWAS) is a study design whereby the effect of a variant, or set of variants is estimated for a large number of traits or outcomes. In TarGene, running a PheWAS is always done implicitely because the outcomes are defined as all variables in the `TRAITS_DATASET` that do not play a particular role in the `ESTIMANDS_FILE`. This file is a description of the target quantities you wish to estimate (as opposed to how you want to estimate them).

To run a PheWAS for 3 genetic variants we will use the `flat` configuration mode. Let's create an `ESTIMANDS_FILE` called `new_phewas_config.yaml` with the following content.

```yaml
type: flat

estimands:
  - type: ATE

variants:
  - 1:238411180:T:C
  - 3:3502414:T:C
  - 2:14983:G:A

outcome_extra_covariates:
  - "Skin colour"
  - "Cheese intake"
```

As a reminder, the `ESTIMANDS_FILE` describes the study design in TarGene. The `estimands` section is set to include Average Treatment Effects (`ATE`) which correspond to single variant effects (the counterpart to the beta coefficients in a linear model). Because there are 3 variants in the `variants` section, 3 PheWAS will actually be run in paralllel. These will be run for all traits in the `UKB_CONFIG` that are not in the `outcome_extra_covariates` section describing additional covariates used to fit the machine-learning models used by TarGene.

!!! note "Note on variant IDs"
    The variant IDs must match the IDs of your BGEN files. Here they are identified by chr:pos:ref:alt but in your files it may be via the rsID. If this is the case, multi-allelic variants are currently not supported.

We can save the following snippet to fully describe our Nextflow run in `my_phewas.config`:

```conf
params {
    ESTIMANDS_CONFIG = "new_phewas_config.yaml"
    BED_FILES = "test/assets/unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    BGEN_FILES = "test/assets/unphased_bgen/ukb_chr{1,2,3}.{bgen,bgen.bgi,sample}"
    UKB_CONFIG = "test/assets/ukbconfig_small.yaml"
    TRAITS_DATASET = "test/assets/dataset.csv"
    UKB_WITHDRAWAL_LIST = "test/assets/withdrawal_list.txt"
}
```

!!! warning "Genotype Files"
    Note that we need to provide both `BED_FILES` and `BGEN_FILES`.

And the command-line to be run:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.13.0 -profile local -c my_phewas.config -resume
```

When the run is complete, you will find the summary statistics in `results/results.summary.yaml`, this file contains a list of all the estimated genetic effects and will look like:

```yaml
- EFFECT_TYPE: "ATE"
  OUTCOME: J40-J47 Chronic lower respiratory diseases
  TREATMENTS:
    - 1:238411180:T:C: "TT => CT"
    - 1:238411180:T:C: "CT => CC"
  WTMLE_TUNEDXGBOOST_TUNEDXGBOOST:
    PVALUE: 0.534213637582527
    COMPONENTS:
      - PVALUE: 0.7120394552494946
        EFFECT_SIZE: -0.0028105792132926276
      - PVALUE: 0.2898059154380053
        EFFECT_SIZE: 0.011852560327880575
  OSE_TUNEDXGBOOST_TUNEDXGBOOST:
    PVALUE: 0.44589759711939014
    COMPONENTS:
      - PVALUE: 0.2624131730093056
        EFFECT_SIZE: -0.008542591675667614
      - PVALUE: 0.23983511767076413
        EFFECT_SIZE: 0.01318599491911467
```

Since we have performed a PheWAS, the treatment variable of each reported effect should correspond to one of the variants specified in the `ESTIAMNDS_FILE`, here 1:238411180:T:C. The effect of this variant on J40-J47 Chronic lower respiratory diseases, is fully defined by the two allelic changes "TT => CT" and "CT => CC". Each of these changes is (by default) estimated in two ways: using a TMLE estimator (`WTMLE_TUNEDXGBOOST_TUNEDXGBOOST` section) and using a One Step estimator (`OSE_TUNEDXGBOOST_TUNEDXGBOOST` section). For each estimator, a global multidimensional p-value is reported corresponding to the global null hypothesis of no effect. Then for each allelic change, a p-value and effect size are further reported. While is small datasets like this mock dataset those two estimators may give widely different results, in large samples they should be very similar.

For more on TarGene's outputs, see [Understanding TarGene's Outputs](@ref).