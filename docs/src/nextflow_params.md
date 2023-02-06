# Index of the pipeline parameters

Here is a list of all the pipeline parameters:

## UK-Biobank Data

- `ENCRYPTED_DATASET`: Path to a UK-Biobank encrypted main dataset.
- `ENCODING_FILE`: If an encrypted dataset is provided, an encoding file must accompany it.
- `DECRYPTED_DATASET`: Path to a UK-Biobank decrypted main dataset. If set, `ENCRYPTED_DATASET` is ignored.
- `TRAITS_CONFIG`: Configuration file describing which traits should be extracted from the main dataset.
- `WITHDRAWAL_LIST`: List of participants withdrawn from the study.
- `QC_FILE`: Genotyping quality control file from the UK-Biobank study.
- `UKBB_BED_FILES`: Path expression to PLINK BED files.
- `UKBB_BGEN_FILES`: Path expression to iputed BGEN files.

## Statistical Parameters

- `PARAMETER_PLAN`: One of "FROM_PARAM_FILES", "FROM_ACTORS". See the section.

If `PARAMETER_PLAN`="FROM_PARAM_FILES":

- `PARAMETER_FILES`: Path expression to the parameter files.

If `PARAMETER_PLAN`="FROM_ACTORS":

- `BQTLS`: CSV file containing binding quantitative trait loci.
- `TRANS_ACTORS`: CSV file containing quantitative trait loci potentially interacting with the previous bqtls.
- `EXTRA_CONFOUNDERS`: Path to additional confounders file, one per line, no header.
- `EXTRA_COVARIATES`: Path to additional covariates file, one per line, no header.
- `ENVIRONMENTALS`: Path to additional environmental treatments file, one per line, no header.
- `ORDERS`: Comma separated list describing the order of desired interaction parameters, 1 for the ATE (no interaction), 2 for pairwise interactions etc... e.g. "1,2"

## TMLE

- `ESTIMATORFILE`: YAML configuration file describing the nuisance parameters learners.
- `PHENOTYPES_BATCH_SIZE`: For a given causal model, phenotypes can be batched to save computational time. Setting this value to 0 will result in max size batches.

## Sieve Variance Correction

- `GRM_NSPLITS`: To fasten GRM computation, it is typically split in batches.
- `NB_VAR_ESTIMATORS`: Number of sieve variance estimates per curve. Setting this value to 0 results in skipping sieve variance correction.
- `MAX_TAU`: Variance estimates are computed for tau ranging from 0 to MAX_TAU
- `PVAL_SIEVE`: To save computation time and disk, only parameters with a p-value below this threshold are considered for sieve variance correction.

## Miscellaneous

- `CALL_THRESHOLD`: Variants from BGEN files are called beyond that threhsold.
- `POSITIVITY_CONSTRAINT`: Treatment variables rarest configuration should have at least that frequency.
- `MAF_THRESHOLD`: Only variants with that minor allele frequency are considered
- `NB_PCS`: Number of principal components to include in confouding adjustment
- `FLASHPCA_EXCLUSION_REGIONS`: Special file from [flashpca](https://github.com/gabraham/flashpca) repository.
- `OUTDIR`: Output directory
