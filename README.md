# TL-Pipeline

> **Warning**:
>   This is rapidly evolving software that should be considered unstable.

Here we provide the main workflow for estimating genetic variants interactions responsible of traits in the UK Biobank using the Targeted Learning framework. For that purpose, we rely on [NextFlow](https://www.nextflow.io/), a software that helps the development of complex parallel and reactive workflows on clouds and clusters.

## Running the Workflow

Please refer to the main [NextFlow](https://www.nextflow.io/) documentation for general usage. The main point being that, depending on your cluster specifications, you will need to provide a specific `myprofile` configuration file. If you are part of the University of Edinburgh and simply using Eddie, then the `eddie` profile is already defined. Then simply run:

```bash
nextflow run TARGENE/targene-pipeline -profile myprofile -resume
```

In order to configure the pipeline to your project you must write a `nextflow.config` file that will provide the parameters described below.

## Configuration

Here is a detailed description of the parameters expected by the pipeline.

### Data source configuration

Currently only the UK-Biobank is supported.

#### UK-Biobank

The UK-Biobank is composed of both genetic data (.bed and .bgen files) and trait data.

- The trait data: The first option is to provide the path to the encrypted dataset `ENCRYPTED_DATASET` which must be decoded via the ukbconv software and the encoding file `ENCODING_FILE`. The second option is to decrypt the dataset only once outside of the pipeline and then use the `DECRYPTED_DATASET` as an input to the pipeline. Finally, since one is usually not interested in all of the traits, and those traits can play a different role in the causal model, the `TRAITS_CONFIG` file described [here](https://github.com/TARGENE/UKBMain.jl) has to be provided.


- The genetic data: We are currently using both .bgen and .bed files. Those are respectively provided with the `UKBB_BGEN_FILES` and `UKBB_BED_FILES` parameters. Since the UK-Biobank genotypic data is split in chromosomes, it should be of the form `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}` and `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}` respectively.

Additional UK-Biobank required files for preprocessing and filtering are:

- `QC_FILE`: A path to the UK-Biobank SNP quaility control `ukb_snp_qc.txt` file.
- `WITHDRAWAL_LIST`: A path to the withdrawal sample list to exclude removed participants from the study.

### Confounding Adjustement configuration

To account for potential confounding effect due to population stratification, we extract principal components from the genetic data using [flashpca](https://github.com/gabraham/flashpca). We follow the recommended procedure for this tool which implies some preprocessing and filtering. The following arguments are compulsory:

- `LD_BLOCKS`: A path to pre-identified linkage desequlibrium blocks around the variants that will be queried for causal effect estimation. Those will be removed from the data.
- `FLASHPCA_EXCLUSION_REGIONS`: A path to the flashpca special exclusion regions which is provided in their repository.
- `NB_PCS` (default: 6): The number of PCA components to extract.

### Parameters configuration files

There are two ways to specify the parameters of interest to your research via the `MODE` parameter. In any case, the goal is to produce a set of parameter files as described in [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl). Those modes are described below:
- `ASBxTransActors`: It is assumed that an initial set of variants has been pre-identified from a previous allele-specific binding (ASB) study as output by the [ball-nf](https://git.ecdf.ed.ac.uk/oalmelid/baal-nf) pipeline: `ASB_FILES` parameter. It is also assumed that another set of potential trans-actors is given in a .csv file: `TRANS_ACTORS_FILE` parameter. In that scenario, the target parameter will be the Interaction Average Treatment Effect between every pair of SNPs. Additionally, if template parameters configuration files containing extra treatments are provided (via `PARAMETER_FILES`), nth-order interaction parameters will be generated.

- `GivenParameters`: This is useful if only a few parameters are of interest and can be written "by hand" or generated outside of the pipeline. Those parameter files can be specified by: `PARAMETER_FILES`.

### Estimators configuration

Almost the last step of the pipeline: targeted estimation. This step uses machine learning estimators which have to be specified via the `ESTIMATORFILE` parameter. This file is more widely described in [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl).

Since the same estimator for `p(T|W)` can be used for multiple target parameters, it may be useful to batch phenotypes using `PHENOTYPES_BATCH_SIZE`(default: 1) in order to reduce the computational burden.


### Sieve Variance Estimation

Finally, the variance estimator can be adjusted via the Sieve Variance Plateau method. For that, we need to compute the GRM which is typically split via `GRM_NSPLITS` (default: 100). Then the number of estimators to compute in the interval [0, `MAX_TAU` (default: 0.8)] is given by `NB_VAR_ESTIMATORS` (default: 0). If `NB_VAR_ESTIMATORS` is set to 0, the Sieve Variance Plateau method will not be applied. It is also possible, in order to reduce the computational burden to perform this correction only if the initial p-value is below a specific threshold `PVAL_SIEVE` (default: 0.05). This is because in the correction will only increase the variance estimate.

