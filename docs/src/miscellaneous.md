# Tweaking additional behaviour

Further Nextflow parameter affecting the behaviour of the pipeline but that does not fit in any previously described category is listed here:

- `CALL_THRESHOLD` (optional, default: 0.9): For putative causal variants (listed in the parameter files described in the [Describing the causal parameters of interest](@ref) section). If a individual's allele's probability is greater than the threshold, then it is called, otherwise it is considered missing.
- `BATCH_SIZE` (optional, default: 400): The set of parameters to be estimated is batched and the TMLE processes will run in parallel across batches on your platform.
- `GENOTYPES_AS_INT` (optional, default: false): If the genotypes should be encoded as a count of the minor allele (0, 1, 2), otherwise the string representation is used.
- `SAVE_IC` (optional, default: true): For all parameters with an p-value below `PVAL_THRESHOLD`, the influence curve is saved. Make sure to keep to `true` if you want to use sieve variance correction, i.e. if `NB_VAR_ESTIMATORS` != 0.
- `OUTDIR` (optional, default: "results"): Output directory
