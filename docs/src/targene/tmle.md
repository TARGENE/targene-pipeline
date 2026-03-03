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

## Addressing Sampling Bias

In population genetics studies, it is possible that the observed population does not adequately represent the population from which it was sampled with respect to the frequency of binary traits like disease. This occurs frequently in case-control studies for which diseased participants (cases) are over-recruited for the study to increase power.

### Why TarGene's Estimands Are Affected

TarGene estimates marginal effects such as the Average Treatment Effect (ATE), which depend on the distribution of covariates in the population. Unlike conditional effects (e.g., conditional odds ratios reported from logistic regression), marginal estimands are not invariant to the sampling strategy. This means that if the covariate distribution in the observed sample differs from the target population, the estimated effect will reflect the sample rather than the true population.

Consider a disease with 1% prevalence in the general population, but a case-control study recruits 50% cases and 50% controls. The ATE estimated from this sample would be computed over a population where half the individuals have the disease-associated covariate profile which is substantially different from the true target population. This can lead to biased estimates of the population-level effect.

### Case-Control Weighted Estimation

To obtain estimates representative of the true population, TarGene supports case-control weighted estimation based on the methodology developed by [Rose and van der Laan (2008)](https://doi.org/10.2202/1557-4679.1115). By providing the known prevalence of the binary phenotype in the target population, we can reweight the empirical distribution to match the population of interest, thereby correcting for sampling bias introduced by the study design.

This procedure involves:
1. **Assigning weights** ``q_0`` to cases and ``(1-q_0)\frac{1}{J}`` to controls, where ``q_0`` is the true population prevalence of the trait of interest and ``J`` is the integer ratio of controls to cases.   
2. **Weighting the nuisace functions**: Both ``Q_Y`` and ``G`` are fitted using case-control weights.
3. **Weighting the targeting step**: The semi-parametric update in (w)TMLE is also appropriately weighted.

A worked example of this bias and the weighted correction procedure can be found [here](https://targene.github.io/TMLE.jl/stable/examples/case_control_experiment/).

!!! warning "Limitations"
    There is no support for case-control weighting of the One-Step Estimator (OSE). Additionally, if a custom machine learning algorithm is used for ``Q_Y`` or ``G``, it must support sample weights (see [MLJ supports weights](https://juliaai.github.io/MLJ.jl/stable/weights/)).
