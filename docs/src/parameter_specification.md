# Parameter plans

## Parameter Files

In this section, by parameter, we mean the statistical parameter that represents the scientific quantity of interest and will be estimated via TarGene. The complete specification of a parameter requires the description of a causal model which can be represented by the following graph.

```@raw html
<div style="text-align:center">
<img src="assets/causal_graph.png" alt="Causal Model" style="width:400px;"/>
</div>
```

Parameters are specified by a YAML file with one section for each variable in the causal model and a section for the parameters that need to be estimated from this causal model. Here is a succinct description for each section and the behaviour of the pipeline for each variable:

- `Treatments`: The treatment variables, typically one or multiple SNPs and potential environmental exposures.
- `Confounders`: Confounding variables, typically the principal components which are computed by the pipeline. This section can (must) be ommited.
- `Covariates`: Additional covariates for the prediction of the traits. Not yet working and must be omitted.
- `Targets`: The traits of interest. The algorithm will loop through all traits in this section for the given treatments, confounders and covariates. This enables the reuse of the propensity score estimation. If this section is omitted (usually the case) all traits will be used.
- `Parameters`: For each target in `Targets` multiple parameters may be of interest depending on the exact case/control scenario of the treatment variables. For instance, since each genotyped locus can take up to 3 different values there could be up to 3 Average treatment Effect parameters if the treatment section consists of only one SNP.

Since an example may be worth a thousand words, here are a couple of such files.

### Average Treatment Effect

One SNP `RSID_10` and 2 ATEs for all traits under study.

```yaml
Treatments:
  - RSID_10
Parameters:
  - name: GG_TO_AG
    RSID_10:
      control: GG
      case: AG
  - name: GG_TO_AA
    RSID_10:
      control: GG
      case: AA
```

### Interaction Average Treatment Effect

Two interacting SNPs `RSID_10` and `RSID_100` with only one target `PHENOTYPE_1` and only one parameter.

```yaml
Treatments:
  - RSID_10
  - RSID_100

Targets:
  - PHENOTYPE_1

Parameters:
  - name: IATE
    RSID_10:
      control: GG
      case: AG
    RSID_100:
      control: GG
      case: AG
```

## Parameter Plans

There are two main ways one can specify parameters that need to be estimated during a targene-pipeline run. This is done via the `MODE` parameter.

### `PARAMETER_PLAN` = `FROM_PARAM_FILES`

This is the most general setting and should match the needs of any project, however it requires some preliminary work. In this setting, one typically provides a set of parameter files as described above. The path to those parameters is then provided with the `PARAMETER_FILES` nextflow parameter. It is not necessary to declare the principal components (from PCA) in the `W` section of those parameter files, they will dynamically be added. The same is true for the target section (`Y`), all variables from the traits dataset that are not part of other sections will be added to the target section. If the `Y` is provided for some parameter file, only those targets will be considered for this specific parameter file.

### `PARAMETER_PLAN` = `FROM_ACTORS`

In this setting the goal is to infer the interaction effect between multiple variants and potential external factors, interacting together via a specific biological mechanism. Typically, multiple sets of variants are of interest and each set is identified with a specific molecule, contributing to the mechanism. In particular, it is assumed that a set of variants, usually binding quantitative trait loci (bQTLs) play a pivotal role because they can be precisely located in the genome and don't suffer from linkage disequilibrium. All interactions of interest are thus defined with respect to that set of genetic variations. Let's Consider the following scenario: we know that a transcription factor binds to molecules `x` and `y` and then differentially binds to specific regions in the genome (`bQTLs`) to regulate downstream genes. We suspect that an alteration of this mechanism is responsible for some diseases. A set of `xQTLs`, associated with the expression of `x` and a set of `yQTLs` associated with the expression of `y` have been identified. Together `xQTLs` and `yQTLs` variants are termed "trans actors". We further suspect that some environmental factors may influence this process. From that scenario, there are many questions that can be asked, for instance : "What is the interaction effects of a bQTL with an environmental factor?". This is a simple pairwise interaction setting and more complex scenarios can be envisaged as described in the following graph.

```@raw html
<div style="text-align:center">
<img src="assets/from_actors.png" alt="FromActors" style="width:800px;"/>
</div>
```

Let us now turn to the pipeline specification for this parameter plan:

- `BQTLS`: A path to a `.csv` file containing at least an `ID` column for each rsID and an optional `CHR` column for the chromosome on which the SNP is located.
- `TRANS_ACTORS`: A path prefix to a set of `.csv` files identifying different trans-acting variants. Each file has the same format as for the `bQTLs`.
- `ENVIRONMENTALS`: A path to a `.txt` file containing a list of environmental exposures with no header and one exposure per line. Each exposure should be available from the trait dataset.
- `ORDERS`: A comma separated string that specifies the various interaction orders of interest. All combinations satisfying the positivity constraint will be generated. The order 1 corresponds to the Average Treatment Effect (ATE) for `bQTLs`, any higher order corresponds to the Interaction Average Treatment Effect (IATE) for the various actors. For example, in the previous scenario, assume we provided `ORDERS`=`1,2`. This would generate parameter files for the estimation of all:
  - ATEs parameters for all bQTLs
  - IATEs parameters for all (bQTLs, xQTLs), (bQTLs, yQTLs), (bQTLs, Envs) pairs.

## Parallelization

Since the same estimator for `p(T|W)` can be used for multiple target parameters, it may be useful to batch phenotypes using `PHENOTYPES_BATCH_SIZE`(default: 1) in order to reduce the computational burden.
