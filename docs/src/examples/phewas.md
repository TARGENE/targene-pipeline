# PheWAS

A phenome-wide association study (PheWAS) is a study design whereby the effect of a variant, or set of variants is estimated for a large number of traits or outcomes. In TarGene, running a PheWAS is always done implicitely because the outcomes are defined as all variables in the `TRAITS_DATASET` that do not play a particular role as per the `ESTIMANDS_FILE`.

To run a PheWAS for 3 genetic variants we will use the `flat` configuration mode. Let's create an `ESTIMANDS_FILE` called "new_phewas_config.yaml" with the following content.

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

The `estimands` section is set to include Average Treatment Effects (`ATE`) which correspond to single variant effects. Because there are 3 variants in the `variants` section, 3 PheWAS will actually be run in paralllel. These will be run for all traits in the `UKB_CONFIG` that are not in the `outcome_extra_covariates` section.

!!! note "Note on variant IDs"
    The variant IDs must match the IDs of your BGEN files. Here they are identified by chr:pos:ref:alt but in your files it may be via the rsID. If this is the case, multi-allelic variants are currently not supported.

The `nextflow.config` file for this study is:

```conf
params {
    ESTIMANDS_CONFIG = "new_phewas_config.yaml"

    // UK-Biobank specific parameters
    BED_FILES = "unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    BGEN_FILES = "unphased_bgen/ukb_chr{1,2,3}.{bgen,bgen.bgi,sample}"
    UKB_CONFIG = "ukbconfig_small.yaml"
    TRAITS_DATASET = "dataset.csv"
    UKB_WITHDRAWAL_LIST = "withdrawal_list.txt"
}
```

!!! warning "Genotype Files"
    Note that we need to provide both `BED_FILES` and `BGEN_FILES`.

And the command-line to be run:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.12.0 -profile local
```
