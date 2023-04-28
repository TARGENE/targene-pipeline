# Index of the pipeline parameters

Here is a list of all the pipeline parameters:

## [Setting a data source](@ref)

- **`ENCRYPTED_DATASET` (required unless a `DECRYPTED_DATASET` is given)**: Path to a UK-Biobank encrypted main dataset.
- **`ENCODING_FILE` (required unless a `DECRYPTED_DATASET` is given)**: If an encrypted dataset is provided, an encoding file must accompany it.
- `DECRYPTED_DATASET` (optional): Path to a UK-Biobank decrypted main dataset. If set, `ENCRYPTED_DATASET` is ignored.
- **`TRAITS_CONFIG` (required)**: Configuration file describing which traits should be extracted from the main dataset.
- **`WITHDRAWAL_LIST` (required)**: List of participants withdrawn from the study.
- **`QC_FILE` (required)**: Genotyping quality control file from the UK-Biobank study.
- **`UKBB_BED_FILES` (required)**: Path expression to PLINK BED files.
- **`UKBB_BGEN_FILES` (required)**: Path expression to iputed BGEN files.

## [Adjusting for confounders](@ref)

- **`LD_BLOCKS` (required)**: A path to pre-identified linkage desequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data.
- **`FLASHPCA_EXCLUSION_REGIONS` (required)**: A path to the flashpca special exclusion regions which is provided in their repository.
- `MAF_THRESHOLD` (optional): Only variants with that minor allele frequency are considered
- `NB_PCS` (optional, default: 6): The number of PCA components to extract.

## [Describing the causal parameters of interest](@ref)

- **`PARAMETER_PLAN` (required, default: "FROM_PARAM_FILES")**: One of "FROM_PARAM_FILES", "FROM_ACTORS".

If `PARAMETER_PLAN`="FROM_PARAM_FILES":

- `PARAMETER_FILES` (required): Path expression to the parameter files.

If `PARAMETER_PLAN`="FROM_ACTORS":

- **`BQTLS` (required)**: A CSV file containing binding quantitative trait loci.
- **`TRANS_ACTORS` (required)**: A prefix to CSV files containing quantitative trait loci potentially interacting with the previous bqtls.
- `EXTRA_CONFOUNDERS` (optional, default: nothing): Path to additional confounders file, one per line, no header.
- `EXTRA_COVARIATES` (optional, default: nothing): Path to additional covariates file, one per line, no header.
- `ENVIRONMENTALS` (optional, default: nothing): Path to additional environmental treatments file, one per line, no header.
- `ORDERS` (optional, default: "1,2"): Comma separated list describing the order of desired interaction parameters, 1 for the ATE (no interaction), 2 for pairwise interactions etc... e.g. "1,2"

## [Specifying a Targeted Estimator](@ref)

- **`ESTIMATORFILE` (required)**: YAML configuration file describing the nuisance parameters learners.
- `POSITIVITY_CONSTRAINT` (optional, default: 0.01): Treatment variables rarest configuration should have at least that frequency.
- `SAVE_IC` (optional, default: true): For all parameters with an p-value below `PVAL_SIEVE`, the influence curve is saved. Make sure to keep to `true` if you want to use sieve variance correction, i.e. if `NB_VAR_ESTIMATORS` != 0.

## [Correcting for population relatedness](@ref)

- `GRM_NSPLITS` (optional, default: 100): To fasten GRM computation, it is typically split in batches.
- `NB_VAR_ESTIMATORS` (optional, default: 0): Number of sieve variance estimates per curve. Setting this value to 0 results in skipping sieve variance correction.
- `MAX_TAU` (optional, default: 0.9): Variance estimates are computed for tau ranging from 0 to MAX_TAU
- `PVAL_SIEVE` (optional, default: 0.05): To save computation time and disk, only parameters with a p-value below this threshold are considered for sieve variance correction.

## [Tweaking additional behaviour](@ref)

- `CALL_THRESHOLD` (optional, default: 0.9): For putative causal variants (listed in the parameter files described in the [Describing the causal parameters of interest](@ref) section). If a individual's allele's probability is greater than the threshold, then it is called, otherwise it is considered missing.
- `BATCH_SIZE` (optional, default: 400): The set of parameters to be estimated is batched and the TMLE processes will run in parallel across batches on your platform.
- `GENOTYPES_AS_INT` (optional, default: false): If the genotypes should be encoded as a count of the minor allele (0, 1, 2), otherwise the string representation is used.
- `OUTDIR` (optional, default: "results"): Output directory
