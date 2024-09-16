# Some Runtime considerations

Targeted Learning can quickly become computationally intensive compared to traditional parametric inference. Here, we describe the drivers of complexity, explain the tricks we use to keep it under control and illustrate with two typical study designs.

## The Drivers of Complexity

Remember that for each estimand of interest, a Targeted Learning estimator requires 3 ingredients:

- An estimator for the propensity score: ``G = P(T|W)``.
- An estimator for the outcome's mean: ``Q_Y = \mathbb{E}[Y|T, W]``.
- A targeting (or debiasing) step.

### Complexity Driven by the Targeting Step

The targeting step is the essential ingredient that makes the estimators in TarGene "work". However, all semi-parametric estimators work slightly differently and are thus not computationally equal. For instance the targeting step of the Targeted Minimum-Loss Estimator (TMLE) requires a generalised linear model fit while the One-Step Estimator (OSE) only requires a simple computation of the bias. The benefit is that, unlike the OSE, the TMLE is guaranteed to always respect the natural bounds of the genetic effect (think probabilities between 0 and 1).

### Complexity Driven by ML Models

The more complex the machine-learning models, the longer the runtime. For instance using a generalized linear model for both ``G`` and ``Q_Y`` would only roughly double the runtime of a regular estimation strategy based on linear models. In contrast, a Super Learning strategy relies on cross-validation and multiple models fits, and is thus much more expensive. However, because the ``G`` model does not involve the outcome ``Y``, it can be efficiently reused hence reducing computational burden if many outcomes are of interest. This is particularly interesting in the PheWAS setting where a single ``G`` model can be used throughout for all estimands.

Similarly, but to a smaller extent, the above targeting step is specific to the variant's genotype change under investigation (e.g., ``CC \rightarrow CT``). Most of the time, multiple genotype changes are of interest (e.g., ``CC \rightarrow CT`` and ``CT \rightarrow TT``) and both ``Q_Y`` and ``G`` can be reused across these changes.

### Complexity Driven by the Estimators

All estimators exist in both a canonical and cross-validated version. The latter relaxes some non-parametric conditions under which the estimators will be asymptotically unbiased. However, it comes at a cost since it requires splitting the data into K folds and fitting the ``G`` and ``Q_Y`` models K times instead of once. Furthermore, the ``G`` model cannot be reused anymore because this outer cross-validation scheme is typically stratified by the outcome as well.

## Controlling Runtime via Parallelisation

If you are running TarGene on a high performance computing platform, you have access to many nodes that can be leveraged. For that purpose we use batching, that is we split all the run's estimands in "cleverly" organised batches to maximize caching of machine learning models. This is controlled by the `BATCH_SIZE (default: 50)` parameter which can be adjusted based on your specific study and platform.

!!! tip "Estimation Resources"
    The estimation is likely going to be the rate limiting step of your run. Each estimation process can further be adjusted in your config file based on your specific problem. For example as follows to retry with more memory as follows. 
    ```config
    process{
        withName: TMLE {
            memory = { 10.GB * task.attempt  }
            time = { 48.hour }
            cpus = 1
        }
    }
    ```

## Examples Through Study Designs

The numbers provided below were obtained with TarGene v0.9.0 and better performance is expected in the future.

### The PheWAS study design

In a PheWAS, one is interested in the effect of a genetic variation across many outcomes (typically around 1000). Because the treatment variable is always the same, the propensity score ``G`` can be reused across all parameters, which drastically reduces computational complexity.

With this setup in mind, the computational complexity is mostly driven by the specification of the learning algorithms for ``Q_Y``, which will have to be fitted for each outcome. In the table below are presented some runtimes for various specifications of ``G`` and ``Q_Y`` using a single cpu. The "Unit runtime" is the average runtime across all estimands averaged over 10 traits and can roughly be extrapolated to bigger studies.

| Estimator | Unit runtime (s) | Extrapolated runtime to 1000 outcomes |
| --- | :---: | :---: |
| `glm` | 4.65 | ≈ 1h20 |
| `glmnet` | 7.19 | ≈ 2h |
| `G-superlearning-Q-glmnet` | 50.05| ≈ 13h45 |
| `superlearning` | 168.98 | ≈ 46h |

Depending on the available resources, this means one can probably afford to use more expensive ML models. This is because the above does not leverage any sort of parallelisation.

### The GWAS study design

In a GWAS, the outcome variable is held fixed and we are interested in the effects of very many genetic variations on this outcome (typically 800 000 for a genotyping array). The propensity score cannot be reused across parameters resulting in a more expensive run.

In this case we look at 3 different genetic variations and only one outcome. In the table below are presented some runtimes for various specifications of ``G`` and ``Q_Y`` using a single cpu. The "Unit runtime" is the average runtime across all estimands and can roughly be extrapolated to bigger studies.

| Estimator file | Continuous outcome unit runtime (s) | Binary outcome unit runtime (s) | Projected Time on HPC (200 folds //) |
| --- | :---: | :---: | :---: |
| `glm` | 5.64 | 6.14 | ≈ 6h30 |
| `glmnet` | 17.46 | 22.24 | ≈ 22h |
| `G-superlearning-Q-glmnet` | 430.54 | 438.67 | ≈ 20 days |
| `superlearning` | 511.26 | 567.72 | ≈ 24 days |

We can see that modern high performance computing platforms definitely enable this study design when using GLMs or GLMNets. It is unlikely however, that you will be able to use Super Learning for any of ``G`` or ``Q_Y`` if you don't have privileged access to such platforms.
