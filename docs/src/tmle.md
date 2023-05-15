# Specifying a Targeted Estimator

The estimator configuration file describes the TMLE specification for the estimation of the parameters defined in the previous section. This YAML configuration file is provided to the pipeline via the `ESTIMATORFILE` parameter.
A set of predifined estimator files are available from [https://github.com/TARGENE/targene-pipeline/estimators/](https://github.com/TARGENE/targene-pipeline/estimators).

That being said, you may need to use a custom set of learning algorithms for your project and write your own and contains 4 sections:

- The `threshold` section is a simple floating point number that specifies the minimum allowed value for p(Treatments|Confounders).
- The `Q_binary`, `Q_continuous` and `G` sections describe the learners for the nuisance parameters. Each of them contains a `model` that corresponds to a valid MLJ model constructor and further keyword hyperparameters. For instance, a Stack can be provided a `measures` argument to evaluate internal algorithms during cross validation. It can also be provided a potentially adaptive `resampling` strategy and the library of `models`. Each of those models can specify a grid of hyperparameters that will individually define a learning algorithm.
  - `Q_binary` corresponds to E[Target| Confounders, Covariates] when the targets are binary variables.
  - `Q_continuous` corresponds to E[Target| Confounders, Covariates] when the targets are continuous variables.
  - `G` corresponds to the joint distribution p(Treatments|Confounders).

Here are two example estimator configurations that can serve as a template.

## Example 1

In this example, Super Learners are used for both `Q` and `G` models. To perform a grid search across model hyperparameters, one can use a list of hyper-parameters. For instance, for the following `Q_continuous` learner, two EvoTreeRegressors will be part of the Super Learner with respectively 10 and 20 trees. The cross-validation procedure can be made adaptive based on the outcome class balance by using the `adaptive: true` option.

```yaml
threshold: 1e-8
# For the estimation of E[Y|W, T]: continuous target
Q_continuous:
  model: Stack
  # Description of the resampling strategy
  resampling:
    type: CV
    # The number of folds is determined based on the data
    adaptive: true
  # List all models and hyperparameters
  models: 
    - type: GLMNetRegressor

    - type: EvoTreeRegressor
      nrounds: [10, 20]

    - type: ConstantRegressor

    - type: HALRegressor
      max_degree: 1
      smoothness_orders: 1
      num_knots: [[10, 5]]
      lambda: 10
      cv_select: false

# For the estimation of E[Y|W, T]: binary target
Q_binary:
  model: Stack
  # Description of the resampling strategy
  resampling:
    type: "StratifiedCV"
    # The number of folds is determined based on the data
    adaptive: true
  # List all models and hyperparameters
  models:
    - type: GLMNetClassifier

    - type: ConstantClassifier

    - type: HALClassifier
      max_degree: 1
      smoothness_orders: 1
      num_knots: [[10, 5]]
      lambda: 10
      cv_select: false

    - type: EvoTreeClassifier
      nrounds: 10

# For the estimation of p(T| W)
G:
  model: Stack
  # Description of the resampling strategy
  resampling:
    type: "StratifiedCV"
    # The number of folds is determined based on the data
    adaptive: true
    # List all models and hyperparameters
  models:
    - type: LogisticClassifier
      fit_intercept: true
    - type: ConstantClassifier
    - type: EvoTreeClassifier
      nrounds: 10
```

## Example 2

While Super Learning is encouraged it is not strictly necessary, here for instance the propensity score is a simple gradient boosting tree and the outcome mode is a logistic regression when the outcome is binary.

```yaml
threshold: 0.001
# For the estimation of E[Y|W, T]: continuous target
Q_continuous:
  model: Stack
  # Description of the resampling strategy
  resampling:
    type: CV
    # The number of folds is determined based on the data
    adaptive: true
  # List all models and hyperparameters
  models: 
    - type: InteractionGLMNetRegressor

    - type: EvoTreeRegressor
      nrounds: [10, 20]
      λ: [0., 1.]
      γ: [0.3]

    - type: ConstantRegressor

# For the estimation of E[Y|W, T]: binary target
Q_binary:
  model: LogisticClassifier
  lambda: 10.

# For the estimation of p(T| W)
G:
  model: EvoTreeClassifier
  nrounds: 10
```

## List of supported models

Here is a list of currently supported models, get in touch if you need more:

| Model       | Regression  | Classification  | Source Package  | Comment |
| ----------- | ----------- | --------------- | --------------- | ------- |
| Gradient Boosting Trees | EvoTreeRegressor | EvoTreeClassifier | [EvoTrees.jl](https://github.com/Evovest/EvoTrees.jl) | Pure Julia implementation of histogram based gradient boosting trees|
| Self Tuning EvoTree | GridSearchEvoTreeRegressor | GridSearchEvoTreeClassifier | [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) | Performs a grid search cross-validation over specified hyper parameters|
| XGBoost | XGBoostRegressor | XGBoostClassifier | [XGBoost](https://xgboost.readthedocs.io/en/stable/) | Julia wrapper around the original libxgboost |
| Self Tuning XGBoost | GridSearchXGBoostRegressor | GridSearchXGBoostClassifier | [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) | Performs a grid search cross-validation over specified hyper parameters|
| Linear Models | LinearRegressor | LogisticClassifier | [MLJLinearModels.jl](https://github.com/JuliaAI/MLJLinearModels.jl) | More models available, see: [the docs](https://juliaai.github.io/MLJLinearModels.jl/stable/api/#MLJ-Interface-1) |
| GLM with penalization constant tuning | GLMNetRegressor | GLMNetClassifier | [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) | This is a simple MLJ API around [GLMNet.jl](https://github.com/JuliaStats/GLMNet.jl), number of cross-validation folds is specified via `nfolds`, parallel execution via `parallel` |
| Same as above with interaction terms | InteractionGLMNetRegressor | InteractionGLMNetClassifier | [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) | order of interaction is specified with `order` |
| Same as above but restricts higher order interactions to specific columns | RestrictedInteractionGLMNetRegressor | RestrictedInteractionGLMNetClassifier | [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl) | order of interaction is specified with `order`, columns are specified by `primary_columns` (a list of columns) and `primary_patterns` (a list of Regular expressions) |
| Highly Adaptive Lasso | HALRegressor | HALClassifier | [HighlyAdaptiveLasso.jl](https://github.com/olivierlabayle/HighlyAdaptiveLasso.jl) | Simple wrapper around the original [R package](https://github.com/tlverse/hal9001) |
| Constant model | ConstantRegressor | ConstantClassifier | [MLJModels.jl](https://github.com/JuliaAI/MLJModels.jl) | Always outputs the target's mean |
