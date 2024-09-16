# Final Tweaks

Further Nextflow parameter affecting the behaviour of the pipeline but that does not fit in any previously described category is listed here:

- `BATCH_SIZE` (optional, default: 50): The set of estimands to be estimated is batched and the TMLE processes will run in parallel across batches on your platform.
- `KEEP_IC` (optional, default: false): For all estimands with a p-value below `PVAL_THRESHOLD`, the influence curve is saved.
- `OUTDIR` (optional, default: "results"): Output directory
