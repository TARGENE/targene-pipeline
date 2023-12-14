# Index of the pipeline parameters

Here is a list of all the pipeline parameters:

## [Setting a data source](@ref)

- **`COHORT` (required):** Current default for this is UKBB. If set to a value other than UKBB, this will not run UKBB-specific trait extraction. If `COHORT` is not UKBB, you must specify your trait data in the `DECRYPTED_DATASET` parameter.
- **`ENCRYPTED_DATASET` (required unless a `DECRYPTED_DATASET` is given)**: Path to a UK-Biobank encrypted main dataset.
- **`ENCODING_FILE` (required unless a `DECRYPTED_DATASET` is given)**: If an encrypted dataset is provided, an encoding file must accompany it.
- `DECRYPTED_DATASET` (optional): Path to a UK-Biobank decrypted main dataset. If set, `ENCRYPTED_DATASET` is ignored. **If you are running this for a non-UKBB cohort, your sample IDs must be specified in the first column of this CSV file, with the column name `SAMPLE_ID`.**
- **`TRAITS_CONFIG` (required)**: Configuration file describing which traits should be extracted from the main dataset.
- **`WITHDRAWAL_LIST` (required)**: List of participants withdrawn from the study.
- `QC_FILE` (optional): Genotyping quality control file from the UK-Biobank study. This is **required** when `COHORT`="UKBB".
- **`BED_FILES` (required)**: Path expression to PLINK BED files.
- **`BGEN_FILES` (required)**: Path expression to imputed BGEN files.

## [Adjusting for confounders](@ref)

- **`LD_BLOCKS` (optional)**: A path to pre-identified linkage disequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data. It is good practice to specify `LD_BLOCKS`, as it will remove SNPs correlated with your variants-of-interest before running PCA. 
- **`FLASHPCA_EXCLUSION_REGIONS` (required)**: A path to the flashpca special exclusion regions which is provided in their repository.
- `MAF_THRESHOLD` (optional): Only variants with that minor allele frequency are considered
- `NB_PCS` (optional, default: 6): The number of PCA components to extract.

## [Study Designs](@ref)

- **`STUDY_DESIGN`** (required, default: "FROM\_PARAM\_FILE"): One of "FROM\_PARAM\_FILE", "FROM\_ACTORS".

If `STUDY_DESIGN`="FROM\_PARAM\_FILE":

- `ESTIMANDS_FILE` (required): Path expression to the parameter files.

If `STUDY_DESIGN`="FROM_ACTORS":

- **`BQTLS` (required)**: A CSV file containing binding quantitative trait loci (bQTLs). If multiple transcription factors (TFs) are included in a single run, you must include a column called `TF`, which specifies the TF associated with each bQTL.
- **`TRANS_ACTORS` (required)**: A prefix to CSV files containing quantitative trait loci potentially interacting with the previous bqtls. If multiple transcription factors (TFs) are included in a single run, you must include a column called `TF`, which specifies the TF associated with each transactor.
- `EXTRA_CONFOUNDERS` (optional, default: nothing): Path to additional confounders file, one per line, no header.
- `EXTRA_COVARIATES` (optional, default: nothing): Path to additional covariates file, one per line, no header.
- `ENVIRONMENTALS` (optional, default: nothing): Path to additional environmental treatments file, one per line, no header.
- `ORDERS` (optional, default: "1,2"): Comma separated list describing the order of desired interaction parameters, 1 for the ATE (no interaction), 2 for pairwise interactions etc... e.g. "1,2"

## [Specifying a Targeted Estimator](@ref)

- **`ESTIMATOR_FILE` (required)**: Julia configuration file describing the nuisance parameters learners.
- `POSITIVITY_CONSTRAINT` (optional, default: 0.01): Treatment variables rarest configuration should have at least that frequency.
- `SAVE_IC` (optional, default: true): For all parameters with an p-value below `PVAL_THRESHOLD`, the influence curve is saved. Make sure to keep to `true` if you want to use sieve variance correction, i.e. if `NB_SVP_ESTIMATORS` != 0.

## [Correcting for population relatedness](@ref)

- `GRM_NSPLITS` (optional, default: 100): To fasten GRM computation, it is typically split in batches.
- `NB_SVP_ESTIMATORS` (optional, default: 0): Number of sieve variance estimates per curve. Setting this value to 0 results in skipping sieve variance correction.
- `MAX_SVP_THRESHOLD` (optional, default: 0.9): Variance estimates are computed for tau ranging from 0 to MAX_SVP_THRESHOLD
- `PVAL_THRESHOLD` (optional, default: 0.05): To save computation time and disk, only parameters with a p-value below this threshold are considered for sieve variance correction.

## [Tweaking additional behaviour](@ref)

- `CALL_THRESHOLD` (optional, default: 0.9): For putative causal variants (listed in the parameter files described in the [Study Designs](@ref) section). If a individual's allele's probability is greater than the threshold, then it is called, otherwise it is considered missing.
- `BATCH_SIZE` (optional, default: 400): The set of parameters to be estimated is batched and the TMLE processes will run in parallel across batches on your platform.
- `OUTDIR` (optional, default: "results"): Output directory
- `RNG` (optional, default: 123): General random seed used where appropriate.

## [Running negative control checks](@ref)

- `MAX_PERMUTATION_TESTS` (optional, default: null): Arbitrarily limits the number of permutation tests performed.
- `PVAL_COL` (optional, default: TMLE_PVALUE): In the output `summary.csv`, the p-value column used to define significant hits.
- `PERMUTATION_ORDERS` (optional, default: "1"): A comma separating string defining the permutation test orders to be performed. (see [Permutation tests](@ref))
- `MAF_MATCHING_RELTOL`(optional, default:0.05): Random variants are chosen with a similar MAF matched with the given relative tolerance.
- `N_RANDOM_VARIANTS`(optional, default: 10): For each hit, that many variants are randomly picked.
