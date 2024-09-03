# Specifying a Targeted Estimator

TarGene is a flexible procedure that does not impose any constraint on the functional form of the relationship between genetic variants, environmental variables and traits. In practice, we rely on [MLJ](https://alan-turing-institute.github.io/MLJ.jl/dev/) for all machine learning algorithms. In population genetics studies, there are two learning algorithms we need to specify:

- ``Q_Y = \mathbb{E}[Y|T, W, C]``: The mean outcome given the treatment, confounders and extra covariates.
- ``G = P(T|W)``: The propensity score, which enables the targeting step of the estimation procedure.

In TarGene, the default it to use for both models, a cross-validated version of [XGBoost](https://xgboost.readthedocs.io/en/stable/) over 3-folds and a range of hyper-parameters.

Furthermore, we estimate genetic effects using both weighted Targeted Minimum-Loss Estimation (wTMLE) and One-Step Estimation (OSE). Note that these estimators are not combined in a single estimator, two distinct estimation procedures are performed. The reason for this is that the finite sample behavior of semi-parametric estimators in population genetics is still under study. Providing the results for both estimation methods thus enables a straightforward comparison.

!!! info "Note"
    Rest assured though, from our experience, these results should be almost indistinguishable. Furthermore, the computational cost of the OSE is almost negligible if (w)TMLE is performed and vice versa.

A default TarGene run will thus result in two estimates for each estimand:

- wTMLE using XGBoost for both ``Q_Y`` and ``G``
- OSE using XGBoost for both ``Q_Y`` and ``G``

Note however that this is not compulsory. A single estimator can be used and the machine-learning models customised. This is done via the `ESTIMATORS_CONFIG` which is either a simple configuration string or a more elaborate Julia file. Both options are further discussed [here](https://targene.github.io/TMLECLI.jl/stable/tmle_estimation/#Specifying-Estimators).

For example, using a configuration string, a computationally cheaper run with `ESTIMATORS_CONFIG=ose--glmnet` uses only One-Step Estimation with a [GLMNet](https://www.jstatsoft.org/article/view/v033i01) model for both ``Q_Y`` and ``G``.
