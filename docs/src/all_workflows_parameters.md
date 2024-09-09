# Index of Workflows Parameters

List of all workflows parameters. For more information on them please visit the dedicated user guides.

## Required Parameters

These must be provided, or the workflow will not run.

- `ESTIMANDS_CONFIG`: YAML configuration file describing the effect sizes of interest.
- `TRAITS_DATASET`: Path to a traits dataset. If you are running this for a non-UKB cohort, your sample IDs must be specified in the first column of this CSV file, with the column name `SAMPLE_ID`.
- `BED_FILES`: Path to PLINK BED files.
- `BGEN_FILES`: Path to **indexed** imputed BGEN files (optional for a GWAS).

## Main Options

These are optional but important as they can have a significant impact on the workflow (speed, estimates, ...).

- `ESTIMATORS_CONFIG (default: wtmle-ose--tunedxgboost)`: Estimator name or Julia file containing the description of the Targeted Estimators to use. To be consistent it should match the argument provided to the previous TarGene run.
- `BATCH_SIZE (default: 50)`: The set of estimands to be estimated is batched and the Targeted Learning processes will run in parallel across batches. This is the main driver of computational speed on High Performance Computing Platforms.
- `COHORT (default: UKB)`: Current default for this is UKB. If set to a value other than UKB, this will not run UKB-specific trait extraction.
- `POSITIVITY_CONSTRAINT (default: 0.01)`: When the list of estimands is generated or validated. Treatment variables' rarest configuration should have at least that frequency. For example if the treatment variables are two variants with minor allele A and T respectively. The rarest configuration will be (AA, TT) and should have a frequency of at least `POSITIVITY_CONSTRAINT`.
- `NB_PCS (default: 6)`: The number of PCA components to extract.

### UK Biobank Specific

For UK Biobank analyses (`COHORT=UKB`).

- `UKB_CONFIG (default: ${projectDir}/assets/ukbconfig.yaml)`: YAML configuration file describing which traits should be extracted and how the population should be subsetted.
- `UKB_ENCODING_FILE`: If the `TRAITS_DATASET` is encrypted, an encoding file must be provided.
- `UKB_WITHDRAWAL_LIST`: List of participants withdrawn from the study.
- `QC_FILE`: Genotyping quality control file from the UK-Biobank study.

## Secondary Options

These are of less importance.

- `SVP (default: false)`: Whether Sieve Variance Plateau correction should be performed.
- `FLASHPCA_EXCLUSION_REGIONS (default: assets/exclusion_regions_hg19.txt)`: A path to the flashpca special exclusion regions.
- `MAF_THRESHOLD (default: 0.01)`: Only variants with that minor allele frequency are considered for PCA.
- `LD_BLOCKS`: A path to pre-identified linkage disequilibrium blocks to be removed from the BED files for PCA. It is good practice to specify `LD_BLOCKS`, as it will remove SNPs correlated with your variants-of-interest before running PCA.
- `VERBOSITY (default: 0)`: Verbosity level of the the Workflow's processes.
- `TL_SAVE_EVERY (default: BATCH_SIZE)`: During the estimation process, results are appended to the file in chunks to free memory.
- `KEEP_IC (default: SVP)`: To save the Influence Curves for each estimate. Depending on the size of your dataset, this can result in very large disk usage.

### Sieve Variance Plateau Options

These are only relevant if you which to apply the experimental sieve variane plateau correction (`SVP=true`).

- **`GRM_NSPLITS`, default: 100**: To fasten GRM computation, it is typically split in batches.
- **`NB_SVP_ESTIMATORS`, default: 100**: Number of sieve variance estimates per curve. Setting this value to 0 results in skipping sieve variance correction.
- **`MAX_SVP_THRESHOLD`, default: 0.9**: Variance estimates are computed for tau ranging from 0 to MAX_SVP_THRESHOLD
- **`PVAL_THRESHOLD`, default: 0.05**: Only results with a p-value below this threshold are considered for Sieve Plateau Variance correction.
- **`ESTIMATOR_KEY`, default: 1**: Identifies an estimator from `ESTIMATORS_CONFIG`. The p-value for `PVAL_THRESHOLD` is computed using the result from this estimator.

## Simulations Specific

These apply to both null and realistic simulations (`-entry NULL_SIMULATION` or `-entry REALISTIC_SIMULATION`).

- `N_REPEATS (default: 2)`: Number of bootstrap samples within a single process for each simulation task.
- `RNGS (default: 1..250)`: Random number generators defining the number of parallel bootstrap processes, each running `N_REPEATS` bootstrap samples.
- `SAMPLE_SIZES (default: [500000])`: The dataset sample sizes for which the simulations will be run.
- `MIN_FACTOR_LEVEL_OCCURENCES (default: 10)`: Each level of each factor should be sampled at least this amount of time.
- `MAX_SAMPLING_ATTEMPTS (default: 10000)`: Number of sampling attempts to be made to respect `MIN_FACTOR_LEVEL_OCCURENCES`.
- `NSAMPLES_FOR_TRUTH (default: 1000000)`: Monte Carlo samples used to evaluate the true effects.

### Realistic Simulation Specific

These are specific to the realistic simulation workflow (`-entry REALISTIC_SIMULATION`).

- `TRAIN_RATIO (default: 6)`: Part of the data (out of 10) which is used for the training set of density estimators.
- `SAMPLE_GA_HITS (default: true)`: Whether additional variants should be sampled from the geneATLAS.
- `GA_MAX_VARIANTS (default: 50)`: The maximum number of variants to sample if `SAMPLE_GA_HITS` is true.
- `GA_DISTANCE_THRESHOLD (default, 1000000)`: Minimum number of base pairs between sampled geneATLAS variants.
- `GA_PVAL_THRESHOLD (default: 1e-5)`: Minimum p-value between geneATLAS sampled variants and the traits.
- `GA_MAF_THRESHOLD (default: 0.01)`: Minimum minor allele frequency of geneATLAS sampled variants.
- `GA_TRAIT_TABLE (default: "${projectDir}/assets/Traits_Table_GeneATLAS.csv")`: Path to the downloaded geneATLAS association table.
