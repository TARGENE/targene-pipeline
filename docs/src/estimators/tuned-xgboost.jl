xgboost_regressor = XGBoostRegressor(tree_method="hist")
xgboost_classifier = XGBoostClassifier(tree_method="hist")

tmle_spec = (
  # Controls caching of data by MLJ machines: turning to `true` may result in faster execution but higher memory usage
  cache        = false,
  # Controls whether the fluctuation is weighted or not
  weighted_fluctuation = false,
  # Propensity score threshold
  threshold    = 1e-8,
  # For the estimation of E[Y|W, T]: continuous target
  Q_continuous = TunedModel(
    model = xgboost_regressor,
    resampling = CV(nfolds=3),
    tuning = Grid(goal=20),
    range = [
        range(xgboost_regressor, :max_depth, lower=3, upper=7), 
        range(xgboost_regressor, :lambda, lower=1e-5, upper=10, scale=:log)
        ],
    measure = log_loss,
    cache=false
    ),
  # For the estimation of E[Y|W, T]: binary target
  Q_binary = TunedModel(
    model = xgboost_classifier,
    resampling = StratifiedCV(nfolds=3),
    tuning = Grid(goal=20),
    range = [
        range(xgboost_classifier, :max_depth, lower=3, upper=7), 
        range(xgboost_classifier, :lambda, lower=1e-5, upper=10, scale=:log)
        ],
    measure = log_loss,
    cache=false
    ),
  # For the estimation of p(T| W)
  G = TunedModel(
    model = xgboost_classifier,
    resampling = StratifiedCV(nfolds=3),
    tuning = Grid(goal=20),
    range = [
        range(xgboost_classifier, :max_depth, lower=3, upper=7), 
        range(xgboost_classifier, :lambda, lower=1e-5, upper=10, scale=:log)
        ],
    measure = log_loss,
    cache=false
)
)