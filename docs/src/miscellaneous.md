# Tweaking additional behaviour

Further Nextflow parameter affecting the behaviour of the pipeline but that does not fit in any previously described category is listed here:

- `CALL_THRESHOLD` (optional, default: 0.9): For putative causal variants (listed in the parameter files described in the [Describing the causal parameters of interest](@ref) section). If a individual's allele's probability is greater than the threshold, then it is called, otherwise it is considered missing.
- `PHENOTYPES_BATCH_SIZE` (optional, default: 0): Depending on the size of your study and constraints of your HPC platform, it may be advantageous to parallelize the TMLE processes further. This can be done by batching phenotypes under study together. Setting this value to 0 will result in max size batches (i.e. no batching).
- `GENOTYPES_AS_INT` (optional, default: false): If the genotypes should be encoded as a count of the minor allele (0, 1, 2), otherwise the string representation is used.
- `OUTDIR` (optional, default: "results"): Output directory