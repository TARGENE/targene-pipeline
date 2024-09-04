# Overview

The **asymptotic** performance of semi-parametric estimators is theoretically optimal. However, in practice, asymptotic regimes are not necessarily achieved even for large sample sizes because some events are extremely rare. This is particularly prevalent in population genetics where some genetic variants and traits are found in less than ``1\%`` of individuals. In the absence of finite sample guarantees, simulation studies provide an effective way to validate statistical methods. An ideal simulation would both yield ground-truth values for the causal effects of interest and be representative of real data. In TarGene, we propose two types of simulations to analyse the performance of semi-parametric estimators.

The first type of simulations, called the null simulation, aims at analysing the behaviour of the estimators when the null hypothesis of "no effect" is true. The rationale behind this simulation is that most genetic variants are believed to have no effect, it is thus of particular importance that the type 1 error rate be controlled appropriately.

The second type of simulations, called the realistic simulation, fits the generating process using data-adaptive flexible machine-learning generative models. As such, it retains important features of the true generating process, but also provides ground truth values via Monte-Carlo sampling.

## Bootstrap Size

Both simulations rely on bootstrap sampling, that is, we resample the dataset ``B`` times to obtain the sampling distribution of the estimators. This is computationally expensive since machine learning models need to be fitted for each bootstrap sample. In order to maximise parallelisation opportunities offered by HPCs, instead of specifying `B`, we specify two parameters: `N_REPEATS (default: 2)` and `RNGS (default: 1..250)`. For each rng in `RNGS`, a single process performing `N_REPEATS` bootstrap will be run. Thus, we effectively have ``B = N_REPEATS \times length(RNGS)``. For example the default uses ``B=500`` bootstrap samples, but if we let `N_REPEATS=10` and `RNGS=1..250`, then ``B=2500``.

## Simulation Tasks

We are interested in the performance of estimators for various estimands, but also potentially for various sample sizes. The triple (sample size, estimator, estimand), thus defines a simulation unit or task. In other words, ``B`` bootstrap samples will be run for each task defined by the cartesian product of all sample sizes, estimators and estimands, hence defining the total simulation scope.

### Sample Sizes

For each sample size in `SAMPLE_SIZES (default: [500000])`, the dataset is resampled with replacement according to the processes defined in the next sections. For example if `SAMPLE_SIZES = [1000, 500000]`, the estimators will be evaluated against all estimands for both sample sizes.

## Estimators

Like for a normal discovery run, the set of estimators to be evaluated is defined by the `ESTIMATORS_CONFIG` parameter, which in this case can be a list. For example if `ESTIMATORS_CONFIG = ["wtmle-ose--glm", "wtmle-ose--tunedxgboost"]`, a total of 4 estimators will actually be evaluated:

- wTMLE with a GLM for both ``Q_Y`` and ``G``
- OSE with a GLM for both ``Q_Y`` and ``G``
- wTMLE with a cross-validated XGBoost for both ``Q_Y`` and ``G``
- OSE with a cross-validated XGBoost for both ``Q_Y`` and ``G``

For more information on defining estimators, have a look at the [Specifying a Targeted Estimator](@ref) section.

## Estimands

Finally, the set of estimands is provided via the `ESTIMANDS_CONFIG` which should be a list of estimands like in the [Custom (Advanced)](@ref) section.

Easier ways to describe estimands may be added in the future but note that given the computational complexity of these simulations, a genome-wide simulation is currently unrealistic.

!!! info "Number of Estimands"
    Keep it small to start with, at the moment experiments have been run with around 30 different genetic effects.
