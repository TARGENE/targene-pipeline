# Description of Simulations Outputs

The output of the simulation workflows is a HDF5 file containing a DataFrame which can be loaded within a [Julia](https://julialang.org/) session for instance. Each row of this DataFrame corresponds to a simulation task, that is a triple (estimand, estimator, sample size). Some of the elements within this DataFrame are objects that are defined within the [TMLE.jl](https://targene.github.io/TMLE.jl/stable/) package and it is recommended to use it for downstream analysis. In particular, if it is not imported before loading the results, the objects will have to be reconstructed by the HDF5 library. The DataFrame contains the following columns:

- ESTIMATOR: The task's estimator
- ESTIMAND: The task's estimand
- SAMPLE_SIZE: The task's sample size
- ESTIMATES: The list of bootstrap estimates, each element is a [TMLE.jl](https://targene.github.io/TMLE.jl/stable/) structure.
- N_FAILED: The number of failed bootstrap samples.
- OUTCOME: The estimand's outcome.
- TRUE_EFFECT: The true genetic effect for the estimand.
- MEAN_COVERAGE: The mean coverage of the estimator for this task.
- MEAN_BIAS: The mean bias of the estimator for this task.
- MEAN_VARIANCE: The mean variance of the estimator for this task.
