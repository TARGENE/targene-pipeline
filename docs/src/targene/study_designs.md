# Defining the Estimands of Interest

The `ESTIMANDS_CONFIG` describes the genetic effects that will be estimated within a TarGene run. In particular, it answers the two following questions:

- What variants are of interest?
- What are the quantities of interest: Average Treatment Effects, Epistatic Interactions, Gene by Environment Interactions, ...?

In all cases, the `ESTIMANDS_CONFIG` is a plain YAML file.

While specifying your estimands, it may be useful to keep the following causal model in mind.

![Causal Graph](../assets/causal_graph.png)

Where ``V_1...V_p`` are a set of genetic variants, the ``Y_1...Y_K`` are a set of traits and ``C`` are a set of additional predictors for ``Y`` but not confounding the genetic effects.

## Preliminaries

In what follows, the genetic effects are defined by the following abbreviations:

- ATE: Average Treatment Effect (G)
- AIE: Average Interaction Effect (GxG, GxE, GxGxE, ...)
- CM: Counterfactual Mean (Unlikely to be of interest)

Also, depending on the study design, some of the YAML file sections are the same, instead of repeating them we list them here.

- `type`: Defines the study design (see below).
- `estimands`: A list of the followings:
  - `type`: The type of generated estimands (`CM`, `ATE`, `AIE`)
  - `orders`: With respect to treatment variables if more than two are provided, the combination orders to be generated. For interactions, the order is always greater or equal than 2.
- `variants`: The list of genetic variants of interest, see below for how this can be specified.
- `extra_treatments`: Environmental treatment variables that are added to the treatments combinations.
- `outcome_extra_covariates`: Additional covariates predictive of the outcomes (`C` in the causal model above).
- `extra_confounders`: Additional confounding variables other than Principal Components (`W` in the causal model above).

!!! info "Outcomes"
    Notice the absence of an `outcomes` section. These are defined by omission, that is, any variable within the trait dataset which is not in the `outcome_extra_covariates` or `extra_confounders` sections will be treated as an outcome.

## Genome-Wide Association Study (GWAS)

This is the most popular study design in population genetics. To run a GWAS, simply provide an `ESTIMANDS_CONFIG` as follows:

```yaml
type: gwas

outcome_extra_covariates:
  - Age At Assessment
```

The `outcome_extra_covariates` section lists additional variables that are predictive of the outcomes in the traits dataset.

## Flat Configuration

In this mode, genetic variants are provided as a flat list, for each estimand type and order specified in the `estimands` section, the estimands are generated using combinations.

For example, with the following file:

```yaml
type: flat

estimands:
  - type: AIE
    orders: [2, 3]
  - type: ATE

variants:
  - RSID_17
  - RSID_99
  - RSID_102

extra_treatments:
  - TREAT_1

outcome_extra_covariates:
  - COV_1

extra_confounders:
  - 21003
  - 22001
```

- Average Treatment Effects of `order` 1 (default) for all (RSID_17, RSID_99, RSID_102, TREAT_1) are generated
- Interaction Effects of `order` 2 and 3 for all combinations of (RSID_17, RSID_99, RSID_102, TREAT_1) are generated. Some of these combinations are (RSID_17, RSID_99), (RSID_102, TREAT_1) or (RSID_102, RSID_99, TREAT_1).

### Grouped Configuration

When the number of variants becomes large, estimating all combinations becomes both computationally expensive and reduces power due to the associated multiple testing burden. If the genetic variants are not chosen at random but based on some biological mechanism (e.g. a transcription factor), it is more efficient to group them. Within each group further subgroups can be defined for the roles played by each variant, only combinations across groups are considered.

For example, the following file describes two groups each consisting of two subgroups:

```yaml
type: groups
estimands:
  - type: AIE
    orders: [2, 3]
variants:
  TF1:
    bQTLs:
      - RSID_17
      - RSID_99
    eQTLs:
      - RSID_102
  TF2:
    bQTLs:
      - RSID_17
      - RSID_198
    eQTLs:
      - RSID_2
extra_treatments:
  - TREAT_1
outcome_extra_covariates:
  - COV_1
extra_confounders:
  - 21003
  - 22001
```

For group `TF1`, only (RSID_17, RSID_102) and (RSID_99, RSID_102) are considered and similarly for group `TF2`.

## Custom (Advanced)

Even though the previous sections should cover most common cases, it may be interesting to define more specific effects. Almost any effect defined within [TMLE.jl](https://targene.github.io/TMLE.jl/stable/) can be estimated in TarGene. It is thus recommended to use the package to create these effects, group them in a `TMLE.Configuration` object and save them to a file using the `TMLE.write_yaml` function.
