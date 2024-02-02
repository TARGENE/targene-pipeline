# Configuration

## Main Outputs

The workflow will produce the following main outputs in the output directory (`OUTDIR`, defaults to `results`):

- An HDF5 file containing estimation results (`HDF5_OUTPUT`, default: results.hdf5)
- An optional JSON file containing estimation results (`JSON_OUTPUT`)
- A Quantile-Quantile summary plot: QQ.png

## Arguments

- **`ESTIMATOR_FILE` (required)**: Estimator name or Julia file containing the description of the Targeted Estimators to use. To be consistent it should match the argument provided to the previous TarGene run.

## Main Options

- **`POSITIVITY_CONSTRAINT` (optional, default: 0.01)**: When the list of estimands is generated or validated. Treatment variables' rarest configuration should have at least that frequency. For example if the treatment variables are two variants with minor allele A and T respectively. The rarest configuration will be (AA, TT) and should have a frequency of at least `POSITIVITY_CONSTRAINT`.
- **`PVAL_THRESHOLD` (optional, default: 0.05)**: Only results with a p-value below this threshold are considered for permutation testing.
- **`ESTIMATOR_KEY` (optional, default: TMLE)**: The p-value for `PVAL_THRESHOLD` is computed using the result from this estimator.

## Secondary Options

- **`RNG` (optional, default: 123)**: General random seed used to induce permutation.
- **`VERBOSITY` (optional, default: 0)**: Verbosity level of the the Workflow's processes.
- **`BATCH_SIZE` (optional, default: 400)**: The set of estimands to be estimated is batched and the Targeted Learning processes will run in parallel across batches.
- **`TL_SAVE_EVERY` (optional: default: 50)**: During the estimation process, results are appended to the file in chunks to free memory.
- **`KEEP_IC` (optional)**: To save the Influence Curves for each estimate. Depending on the size of your dataset, this can result in massive disk usage.
