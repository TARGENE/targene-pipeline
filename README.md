# TL-Pipeline

Here we provide the main workflow for estimating genetic variants interactions responsible of traits in the UK Biobank using the Targeted Learning framework. For that purpose, we rely on [NextFlow](https://www.nextflow.io/), a software that helps the development of complex parallel and reactive workflows on clouds and clusters.

## Running the Workflow

Please refer to the main [NextFlow](https://www.nextflow.io/) documentation for general usage. The main point being that, depending on your cluster specifications, you will need to provide a specific `myprofile` configuration file. If you are part of the University of Edinburgh and simply using Eddie, then the `eddie` profile is already defined. Then simply run:

```bash
nextflow run TARGENE/targene-pipeline -profile myprofile -resume
```

In order to configure the pipeline to your project you must write a `nextflow.config` file that will provide the parameters described below.

## Configuration

Here is a detailed description of the parameters expected by the pipeline.

### Data sources

Currently only the UK-Biobank is supported.

#### UK-Biobank

The following arguments are required:

- `UKBMAIN_DATASET`: It is the dataset containing individuals' traits, typically extracted using the `ukbconv` program provided by the UK-Biobank. Traits are extracted and split based on the role they play in the causal model according to a configuration file described [here](https://github.com/TARGENE/UKBMain.jl). The configuration file also enables filtering of individuals.
- `UKBB_BGEN_FILES`: The imputed genotypes from the UK-Biobank. Since the UK-Biobank genotype data is split in chromosomes, it should be of the form `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}`.
- `UKBB_BED_FILES`: The base genotypes in PLINK BED format. Again in the form of `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}`.
- `QC_FILE`: A path to the UK-biobank SNP quaility control `ukb_snp_qc.txt` file.
- `WITHDRAWAL_LIST`: A path to the withdrawal sample list to exclude removed participants from the study.

### Confounding Adjustement

To account for potential confounding effect due to population stratification, we extract principal components from the genetic data using [flashpca](https://github.com/gabraham/flashpca). We follow the recommended procedure for this tool which implies some preprocessing and filtering. The following arguments are compulsory:

- `LD_BLOCKS`: A path to pre-identified linkage desequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data.
- `FLASHPCA_EXCLUSION_REGIONS`: A path to the flashpca special exclusion regions which is provided in their repository.
- `NB_PCS`: The number of PCA components to extract.

### Queries Generation

There are currently two possible procedures which are specified by the `QUERIES_MODE` argument.
- `QUERIES_MODE` = "ASBxTransActors". It is assumed that an initial set of variants has been pre-identified from a previous allele-specific binding (ASB) study as output by the [ball-nf](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf) pipeline. It is also assumed that another set of potential trans-actors if given in a csv file. The cross product of those two sets of variants is taken to generate the queries. For more information on the procedure and expected format of those files see: `julia bin/generate_queries.jl --help`. Two arguments thus have to be passed:
    - `ASB_FILES`: Path to all files output by the ball-nf pipeline.
    - `TRANS_ACTORS_FILE`: Path to the trans-actors files.
    - Additionaly, the `THRESHOLD` argument can be overided to specify the hard-calling threshold to convert probabilities to genotypes values.

- `QUERIES_MODE` = "given". This is useful if only a few query files have to be generated and can be written "by hand". Then the following argument has to be provided:
    - `QUERY_FILES`: Path to the query files.

### ESTIMATION

Almost the last step of the pipeline: targeted estimation. This step uses "SuperLearning" (Stacking) which means a configuration for this learning step has to be provided:

- `ESTIMATORFILE`: This configuration will be used to build super-learners estimators
- `CROSSVAL`: A boolean to indicate weither an additional cross validation procedure is run to evaluate the models in each super learning step. Unfortunately this is currently run aside of TMLE and thus quite expensive. Only use it if you really need those figures.
- `PHENOTYPES_BATCH_SIZE`: Estimation is parallelized over queries. If the number of queries is low it can be advantageous to parallelize over phenotypes too.


### Sieve Variance Estimation

Arguments for the sieve variance correction step:

- `GRM_NSPLITS`: The number of sub grm parts to be computed since the full GRM requires more than 1TB of memory.
- `NB_VAR_ESTIMATORS`: Number of estimators to compute, the interval [0, MAX_TAU], will be split.
- `MAX_TAU`: Maximum distance up to which individuals will be included in the estimation. Previous analysis showed that above 0.8, some side effects seem to occur. The theoretical maximum is 2.
- `PVAL_SIEVE`: Only traits for which the IID pvalue is lower than this threshold will be considered 
for sieve variance correction. This is because in theory the Sieve Variance curve is supposed to be monotically increasing.

