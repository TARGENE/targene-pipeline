# Specifying a Targeted Estimator

TarGene is a flexible procedure that does not impose any constraint on the functional form of the relationship between genetic variants, environmental variables and traits. In practice, we rely on [MLJ](https://alan-turing-institute.github.io/MLJ.jl/dev/) to provide machine learning algorithms. In population genetics studies, there are two learning algorithms we need to specify:

- `E[Y|T, W, C]`: The mean outcome given the treatment, confounders and extra covariates. It is commonly denoted by `Q` in the Targeted Learning literature. In reality, we will need one specification for continuous outcomes and one specification for binary outcomes.
- `p(T|W)`: The propensity score, which enables the targeting step of the estimation procedure. It is commonly denoted by `G` in the Targeted Learning literature.

There are two main ways to define targeted estimators, from a predefined configuration or from a custom file.

## Predefined estimators

A set of predefined estimators is readily available and can be accessed by using the configuration's name as the `ESTIMATOR_FILE` parameter. For instance, we provide the following:

- G-superlearning-Q-glm
- G-superlearning-Q-glmnet
- glm-with-interactions-for-Q
- glm
- glmnet-with-interactions-for-Q
- glmnet
- superlearning-with-interactions-for-Q
- superlearning
- tuned-xgboost

To see what the estimand files for these examples look like, please refer to the [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) package, located here: https://github.com/TARGENE/TargetedEstimation.jl/tree/main/estimators-configs.

All using TMLE as the meta statistical inference method. While these should cover most use cases, it may be useful to define a custom estimation strategy. This can be achieved by writing a small [Julia](https://julialang.org/) file described below.

## Custom estimators from a file

Writing a [Julia](https://julialang.org/) estimators file is the most flexible way to define an estimation strategy for your study.

This file should simply define a [NamedTuple](https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple) called `ESTIMATORS` and listing estimators to be used for inference. For example, the following:

```julia
default_models = TMLE.default_models(
  # For the estimation of E[Y|W, T]: continuous outcome
  Q_continuous = LinearRegressor(),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary = LogisticClassifier(lambda=0.),
  # For the estimation of p(T| W)
  G = LogisticClassifier(lambda=0.)
)

ESTIMATORS = (
  TMLE_weighted   = TMLEE(models=default_models, weighted=true),
  TMLE_unweighted = TMLEE(models=default_models, weighted=false),
  OSE             = OSE(models=default_models)
)
```

defines three estimators:

1. A weighted-fluctuation Targeted Maximum Likelihood Estimator: `TMLE`
2. An unweighted-fluctuation Targeted Maximum Likelihood Estimator: `TMLE`
3. A One-Step Estimator: `OSE`

All estimators will learn the nuisance functions `Q` and `G` with the provided `models` `NamedTuple`:

- `Q_continuous`: A MLJ model used for the estimation of `E[Y|T, W, C]` when the outcome `Y` is continuous.
- `Q_binary`: A MLJ model used for the estimation of `E[Y|T, W, C]` when the outcome `Y` is binary.
- `G`: A MLJ model used for the estimation of `p(T|W)`.

For the list of available models and resampling strategies, checkout the [TargetedEstimation documentation](https://targene.github.io/TargetedEstimation.jl/stable/models/).

For full details, on available estimators and how to specify them, visit the [TMLE.jl documentation](https://targene.github.io/TMLE.jl/stable/).

## Current recommanedations

At present, there are two main use cases for Targene: (1) investigating an average treatment effect (ATE) or (2) an interaction average treatment effect (IATE). Here we present our current recommendation fo estimators for these two types of analyses.

1. ATE

```julia
default_models = TMLE.default_models(
  # For the estimation of E[Y|W, T]: continuous outcome
  Q_continuous = LinearRegressor(),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary = LogisticClassifier(lambda=0.),
  # For the estimation of p(T| W)
  G = LogisticClassifier(lambda=0.)
)

ESTIMATORS = (
  TMLE_weighted   = TMLEE(models=default_models, weighted=true),
  OSE             = OSE(models=default_models)
)
```

2. IATE

```
xgboost_classifier = XGBoostClassifier(tree_method="hist")
resampling = JointStratifiedCV(patterns=[r"<regex-expression-matching-all-snps-in-bgen-files>"], resampling=StratifiedCV(nfolds=3))
 
default_models = TMLE.default_models(
  # For the estimation of E[Y|W, T]: continuous outcome
  Q_continuous = Pipeline(
    RestrictedInteractionTransformer(order=2, primary_variables_patterns=[r"<regex-expression-matching-all-snps-in-bgen-files>"]),
    GLMNetRegressor(resampling=resampling),
    cache = false
  ),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary    = Pipeline(
    RestrictedInteractionTransformer(order=2, primary_variables_patterns=[r"<regex-expression-matching-all-snps-in-bgen-files>"]),
    GLMNetClassifier(resampling=resampling),
    cache = false
  ),
  # For the estimation of p(T| W)
  G = Stack(
    metalearner        = LogisticClassifier(lambda=0., fit_intercept=false),
    resampling         = resampling,
    cache              = false,
    glmnet             = GLMNetClassifier(resampling=resampling),
    lr                 = LogisticClassifier(lambda=0.),
    tuned_xgboost      = TunedModel(
        model = xgboost_classifier,
        resampling = resampling,
        tuning = Grid(goal=20),
        range = [
            range(xgboost_classifier, :max_depth, lower=3, upper=7),
            range(xgboost_classifier, :lambda, lower=1e-5, upper=10, scale=:log)
            ],
        measure = log_loss,
        cache=false
    )
  )
)
 
ESTIMATORS = (
  TMLE = TMLEE(models=default_models, weighted=true),
  OSE  = OSE(models=default_models)
)
```

The `RestrictedInteractionTransformer` investigates all interactions between your genotypes-of-interest and other variables in your model, without computing the interaction between every PC. You must fill in the `<regex-expression-matching-all-snps-in-bgen-files>` with a REGEX expression that matches all SNP IDs in your BGEN files and additional environmental factors included in a run. Since the genotypes are pulled from the BGEN files, the format of these IDs must be checked in advance, as they may follow the structure of rsIDs or they may not. 

Some examples for `<regex-expression-matching-all-snps-in-bgen-files>` include:
1. "^rs[0-9]+"
2. "^chr([1-9]|1[0-9]|2[0-2]):\d+:[GATC]+:[GATC]+"
3. "(^rs[0-9]+|Genetic-Sex)"

A BioBank-scale simulation study is currently being conducted, and these recommendations will be updated in a later version.