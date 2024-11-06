# Interaction Study

In TarGene you can investigate the interacting effects between genetic variants and between genetic variants and the environment up to any order (GxG, GxE, GxGxG, ...). However, the computational complexity associated with interaction testing grows exponentially. Interaction studies are thus very often targeted, for instance towards a specific transcription factor `TF`. For this transcription factor, multiple variants may have been identified because they are believed to play different roles. In TarGene such a scenario can be encoded via the `group` mode of the `ESTIMANDS_CONFIG` file. An example is more informative than a long description. Let's save the following content to a "new_interaction_config.yaml" file.

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

The `estimands` section describes the types of estimands we are interested in, here Average Interaction Effects (`AIE`). We would like to scan both orders 2 and 3.

The `variants` section contains groups of variants for which the interactions will be generated independently. Here there is only one group called `TF`, and it contains two subgroups: `variants_1` and `variants_2`. This means that combinations of `variants_1` and `variants_2` will define the interactions under investigation. Here they are simply (1:238411180:T:C, 2:14983:G:A) and (3:3502414:T:C, 2:14983:G:A). Why did we set the interaction `orders` to be both 2 and 3 since only pairwise interactions can be generated here?

This is because the `extra_treatments` section can list additional environmental variables, thus increasing the potential interaction order. In this case, there is one: Cheese intake. This means that additional pairwise interactions will be generated: (1:238411180:T:C, Cheese intake), (2:14983:G:A, Cheese intake) and (3:3502414:T:C, Cheese intake). But also the following order 3 interactions: (1:238411180:T:C, 2:14983:G:A, Cheese intake) and (3:3502414:T:C, 2:14983:G:A, Cheese intake).

!!! note "Note on variant IDs"
    The variant IDs must match the IDs of your BGEN files. Here they are identified by chr:pos:ref:alt as but in your case it may be via the rsID.

Finally, the `outcome_extra_covariates` describes additional variables predictive of the outcome, but that do not confound the effect of treatment variables (`variants` and `extra_treatments`).

Importantly, all traits variables in the `TRAITS_DATASET` that are not listed as either `outcome_extra_covariates` or `extra_treatments` will be considered as outcome variables.

We will now point to this newly created `ESTIMANDS_CONFIG` in our `nextflow.config` file:

```conf
params {
    ESTIMANDS_CONFIG = "new_interaction_config.yaml"

    // UK-Biobank specific parameters
    BED_FILES = "unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    BGEN_FILES = "unphased_bgen/ukb_chr{1,2,3}.{bgen,bgen.bgi,sample}"
    UKB_CONFIG = "ukbconfig_small.yaml"
    TRAITS_DATASET = "dataset.csv"
    UKB_WITHDRAWAL_LIST = "withdrawal_list.txt"
}
```

As usual, the pipeline can then be run as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.11.1 -profile local
```
