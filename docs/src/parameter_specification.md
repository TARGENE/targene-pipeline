# Describing the causal parameters of interest

## Parameter Files

In this section, by parameter, we mean the statistical parameter that represents the scientific quantity of interest and will be estimated via TarGene. The complete specification of a parameter requires the description of a causal model which can be represented by the following graph.

```@raw html
<div style="text-align:center">
<img src="../assets/causal_graph.png" alt="Causal Model" style="width:600px;"/>
</div>
```

Parameters are specified by a YAML file containing the list of parameters to be estimated. Each parameter is fully determined by four fields:

- `type`: The parameter type, one of:
  - ATE: Average Treatment Effect
  - IATE: Interaction Average Treatment Effect
  - CM: Conditional Mean
- `target`: The trait of interest, as specified in the `TRAITS_CONFIG` file in the `phenotypes/name` field. You can also use the wildcard "*" to signify that you want to estimate this parameter accross all traits in the dataset.
- `treatment`: The treatment variables and associated control/case settings (see example below).
- `confounders`: A list of confounding variables. Note that principal components will be added to that list by default and must not be provided here. You can provide an empty list.
- `covariates`: This is optional and correspond to a list of additional covariates for the prediction of the trait.

Here is an example parameter file:

```yaml

Parameters:
  - type: IATE
    target: "*"
    treatment: (RSID_10 = (control = "AA", case = "AC), Sun_Exposure = (control = 1, case = 0))
    confounders: [Economic_background]
    covariates: [Age]
  - type: ATE
    target: PHENOTYPE_2
    treatment: (RSID_10 = (control = "AA", case = "CC"), RSID_100 = (control = "GC", case = "CC"))
    confounders: []
  - type: CM
    target: PHENOTYPE_2
    treatment: (RSID_10 = "AA",)
    confounders: []
```

Note that variants can be encoded either with an explicit string representation (e.g. "AC") or via an integer (0, 1, 2) representing the number of minor alleles in the genotype. All variants in a parameter file should however respect the same encoding. String encoded variants are later one-hot-encoded during the TMLE step.

## Parameter Plans

At the moment, there are two main ways one can specify parameters that need to be estimated during a targene-pipeline run. This is done via the `MODE` parameter.

### `PARAMETER_PLAN` = `FROM_PARAM_FILE`

This is the most general setting and should match the needs of any project, however it requires some preliminary work. In this setting, one typically provides a set of parameter files as described above. If you are interested in only a few parameters it may be acceptable to write them by hand. Otherwise it is best to generate them using a programming language (for instance using [TMLE.jl](https://targene.github.io/TMLE.jl/stable/)). The path to those parameters is then provided with the `PARAMETER_FILE` nextflow parameter. See the previous section on [Parameter Files](@ref).

### `PARAMETER_PLAN` = `FROM_ACTORS`

In this setting the goal is to infer the interaction effect between multiple variants and potential external factors, interacting together via a specific biological mechanism. Typically, multiple sets of variants are of interest and each set is identified with a specific molecule, contributing to the mechanism. In particular, it is assumed that a set of variants, usually binding quantitative trait loci (bQTLs) play a pivotal role. All interactions of interest are thus defined with respect to that set of genetic variations. Let's Consider the following scenario: we know that a transcription factor binds to molecules `x` and `y` and then differentially binds to specific regions in the genome (`bQTLs`) to regulate downstream genes. We suspect that an alteration of this mechanism is responsible for some diseases. A set of `xQTLs`, associated with the expression of `x` and a set of `yQTLs` associated with the expression of `y` have been identified. Together `xQTLs` and `yQTLs` variants are termed "trans actors". We further suspect that some environmental factors may influence this process. From that scenario, there are many questions that can be asked, for instance : "What is the interaction effects of a bQTL with an environmental factor?". This is a simple pairwise interaction setting and more complex scenarios can be envisaged as described in the following graph.

```@raw html
<div style="text-align:center">
<img src="../assets/from_actors.png" alt="FromActors" style="width:600px;"/>
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

Since the same estimator for `p(T|W)` can be used for multiple target parameters, it may be useful to batch phenotypes using `BATCH_SIZE`(default: 400) in order to reduce the computational burden.
