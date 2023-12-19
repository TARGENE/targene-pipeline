# Index of the pipeline parameters

Here is a list of all the pipeline parameters:

## [Setting a data source](@ref)

- **`COHORT` (required):** Current default for this is UKBB. If set to a value other than UKBB, this will not run UKBB-specific trait extraction. If `COHORT` is not UKBB, you must specify your trait data in the `TRAITS_DATASET` parameter.
- **`TRAITS_DATASET` (required):** Path to a traits dataset. If you are running this for a non-UKBB cohort, your sample IDs must be specified in the first column of this CSV file, with the column name `SAMPLE_ID`.
- **`BED_FILES` (required)**: Path expression to PLINK BED files.
- **`BGEN_FILES` (required)**: Path expression to imputed BGEN files.
- **`UKB_CONFIG` (required for COHORT=UKBB)**: Configuration file describing which traits should be extracted from the main dataset.
- `UKB_ENCODING_FILE` (optional): If the `TRAITS_DATASET` is encrypted, an encoding file must be provided.
- `WITHDRAWAL_LIST` (optional): List of participants withdrawn from the study.
- `QC_FILE` (optional): Genotyping quality control file from the UK-Biobank study.

## [Adjusting for confounders](@ref)

- `LD_BLOCKS` (optional): A path to pre-identified linkage disequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data. It is good practice to specify `LD_BLOCKS`, as it will remove SNPs correlated with your variants-of-interest before running PCA. 
- `FLASHPCA_EXCLUSION_REGIONS` (optional, default: data/exclusion_regions_hg19.txt): A path to the flashpca special exclusion regions.
- `MAF_THRESHOLD` (optional, default: 0.01): Only variants with that minor allele frequency are considered
- `NB_PCS` (optional, default: 6): The number of PCA components to extract.

## [Study Designs](@ref)

- **`STUDY_DESIGN`** (required, default: "FROM\_PARAM\_FILE"): One of "FROM\_PARAM\_FILE", "FROM\_ACTORS".

If `STUDY_DESIGN`=`CUSTOM`:

- `ESTIMANDS_FILE` (required): Path expression to the estimands file.

If `STUDY_DESIGN`=`ALLELE_INDEPENDENT`:

- `ESTIMANDS_FILE` (required): YAML configuration file.

If `STUDY_DESIGN`=`FROM_ACTORS`:

- **`BQTLS` (required)**: A CSV file containing binding quantitative trait loci (bQTLs). If multiple transcription factors (TFs) are included in a single run, you must include a column called `TF`, which specifies the TF associated with each bQTL.
- **`TRANS_ACTORS` (required)**: A prefix to CSV files containing quantitative trait loci potentially interacting with the previous bqtls. If multiple transcription factors (TFs) are included in a single run, you must include a column called `TF`, which specifies the TF associated with each transactor.
- `EXTRA_CONFOUNDERS` (optional, default: nothing): Path to additional confounders file, one per line, no header.
- `EXTRA_COVARIATES` (optional, default: nothing): Path to additional covariates file, one per line, no header.
- `ENVIRONMENTALS` (optional, default: nothing): Path to additional environmental treatments file, one per line, no header.
- `ORDERS` (optional, default: "1,2"): Comma separated list describing the order of desired interaction estimands, 1 for the ATE (no interaction), 2 for pairwise interactions etc... e.g. "1,2"

## [Specifying a Targeted Estimator](@ref)

- **`ESTIMATOR_FILE` (required)**: Julia configuration file describing the nuisance function learning algorithms.
- `POSITIVITY_CONSTRAINT` (optional, default: 0.01): Treatment variables rarest configuration should have at least that frequency.
- `KEEP_IC` (optional, default: false): For all parameters with an p-value below `PVAL_THRESHOLD`, the influence curve is saved.

## [Correcting for population relatedness](@ref)

- `SVP` (optional, default: false): To perform Sieve Variance Correction.
- `GRM_NSPLITS` (optional, default: 100): To fasten GRM computation, it is typically split in batches.
- `NB_SVP_ESTIMATORS` (optional, default: 100): Number of sieve variance estimates per curve. Setting this value to 0 results in skipping sieve variance correction.
- `MAX_SVP_THRESHOLD` (optional, default: 0.9): Variance estimates are computed for tau ranging from 0 to MAX_SVP_THRESHOLD
- `PVAL_THRESHOLD` (optional, default: 0.05): To save computation time and disk, only parameters with a p-value below this threshold are considered for sieve variance correction.
- `SVP_ESTIMATOR_KEY` (optional, default: "TMLE"): The estimator to use for Sieve Variance Plateau correction.

## [Tweaking additional behaviour](@ref)

- `CALL_THRESHOLD` (optional, default: 0.9): For putative causal variants (listed in the parameter files described in the [Study Designs](@ref) section). If a individual's allele's probability is greater than the threshold, then it is called, otherwise it is considered missing.
- `BATCH_SIZE` (optional, default: 400): The set of parameters to be estimated is batched and the TMLE processes will run in parallel across batches on your platform.
- `OUTDIR` (optional, default: "results"): Output directory
- `RNG` (optional, default: 123): General random seed used where appropriate.

## [Running negative control checks](@ref)

- `MAX_PERMUTATION_TESTS` (optional, default: null): Arbitrarily limits the number of permutation tests performed.
- `ESTIMATOR_KEY` (optional, default: "TMLE"): The estimator from `ESTIMATOR_FILE` to use to asses significance.
- `PERMUTATION_ORDERS` (optional, default: "1"): A comma separating string defining the permutation test orders to be performed. (see [Permutation tests](@ref))
- `MAF_MATCHING_RELTOL`(optional, default:0.05): Random variants are chosen with a similar MAF matched with the given relative tolerance.
- `N_RANDOM_VARIANTS`(optional, default: 10): For each hit, that many variants are randomly picked.
