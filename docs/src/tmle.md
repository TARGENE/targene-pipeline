# Targeted Maximum Likelihood Estimation

The estimator configuration file describes the TMLE specification for the estimation of the parameters defined in the previous section. This YAML configuration file is provided to the pipeline via the `ESTIMATORFILE` parameter and contains 4 sections:

- The `threshold` section is a simple floating point number that specifies the minimum allowed value for p(Treatments|Confounders).
- The `Q_binary`, `Q_continuous` and `G` sections describe the learners for the nuisance parameters. Each of them contains a `model` that corresponds to a valid MLJ model constructor and further keyword hyperparameters. For instance, a Stack can be provided a `measures` argument to evaluate internal algorithms during cross validation. It can also be provided a potentially adaptive `resampling` strategy and the library of `models`. Each of those models can specify a grid of hyperparameters that will individually define a learning algorithm.
  - `Q_binary` corresponds to E[Target| Confounders, Covariates] when the targets are binary variables.
  - `Q_continuous` corresponds to E[Target| Confounders, Covariates] when the targets are continuous variables.
  - `G` corresponds to the joint distribution p(Treatments|Confounders).


Since the same estimator for `p(T|W)` can be used for multiple target parameters, it may be useful to batch phenotypes using `PHENOTYPES_BATCH_SIZE`(default: 1) in order to reduce the computational burden.
