tmle_spec = (
  # Controls caching of data by MLJ machines: turning to `true` may result in faster execution but higher memory usage
  cache        = false,
  # Controls whether the fluctuation is weighted or not
  weighted_fluctuation = false,
  # Propensity score threshold
  threshold    = 1e-8,
  # For the estimation of E[Y|W, T]: continuous target
  Q_continuous = Pipeline(
    RestrictedInteractionTransformer(order=2, primary_variables_patterns=[r"^rs[0-9]+"]),
    GLMNetRegressor(resampling=CV(nfolds=3)),
    cache = false
  ),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary    = Pipeline(
    RestrictedInteractionTransformer(order=2, primary_variables_patterns=[r"^rs[0-9]+"]),
    GLMNetClassifier(resampling=StratifiedCV(nfolds=3)),
    cache = false
  ),
  # For the estimation of p(T| W)
  G           = GLMNetClassifier(resampling=StratifiedCV(nfolds=3))
)