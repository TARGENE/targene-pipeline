# Setup (Read First)

A population genetics analysis via TarGene is a 4 steps process where steps 1 and 2 are done only once.

## 1. Installation

Install TarGene's dependencies as per the [Installation](@ref) section. Here I will assume that both Nextflow and Singularity are installed. To make sure this is the case you can for instance try

```bash
nextflow -v
```

and

```bash
singularity --version
```

## 2. Obtaining the Data

Make sure the data is available. Here we will rely on TarGene's test data located [here](https://github.com/TARGENE/targene-pipeline/tree/main/test/assets). To obtain it, the easiest way is maybe to clone the entire project like so 

```bash
git clone https://github.com/TARGENE/targene-pipeline
```

Then open the repository's folder using your favourite editor. All examples hereafter assume that you have a terminal opened in the repository's root folder.

!!! note "This is only to get the data"
    Even though we have cloned the repository, we do not need to look at the code at all!

## 3. Writing a Run Configuration File

Write a run configuration file for your project, this can be decomposed further into 2 sub-steps:

1. Write the part of your configuration file specific to your platform (HPC, laptop, ...). Here we will simply use the `local` profile which runs every process using singularity and your local cpu. However, you could change this to use docker instead, or resort to a scheduler if you are working on an HPC. It might be worth having a look at [this page](https://nextflow.io/docs/latest/config.html) if you have never worked with Nextflow before. All existing TarGene's configuration files are in the `conf` folder.
2. **Write the part of your configuration file specific to your run**. This is the topic of the remainder of the following sections. That is, each example requires you to write this file. Feel free to experiment by changing the workflow parameters (see [Index of Workflows Parameters](@ref)).

## 4. General Input Parameter Description

- `COHORT`: Here we assume that the cohort is the UK-Biobank dataset (default). If you have already extracted your covariates and phenotypes you can just use `COHORT=CUSTOM` and the below `UKB_CONFIG` won't be used.
- `UKB_CONFIG`: In the examples, we assume that the dataset is a raw UK Biobank dataset which is not readily interpretable by machine-learning models. How does TarGene know how to extract covariates and phenotypes from the the UK Biobank main dataset? This is thanks to the `UKB_CONFIG` file, which maps UK Biobank data fields to traits (see [The `UKB_CONFIG` Configuration File](@ref)). In our case, this file (`test/assets/ukbconfig_gwas.yaml`) contains both the outcome of interest (BMI) and the extra predictors we need (Number of vehicles in household and Cheese intake).

```yaml
traits:
  - fields:
      - "21001"
    phenotypes:
      - name: "Body mass index (BMI)"
  - fields:
      - "728"
    phenotypes:
      - name: "Number of vehicles in household"
  - fields:
      - "1408"
    phenotypes:
      - name: "Cheese intake"
```

Note that the outcomes are defined implicitely, any trait in the `UKB_CONFIG` file which are not in the `outcome_extra_covariates` of the `ESTIMANDS_FILE` (see later) will be considered as outcomes. So you can run multiple GWAS at once by simply adding another trait definition to the above file.

- `TRAITS_DATASET`: This is a CSV file containing both covariates and phenotypes of interest. It should also contain an `eid` column to uniquely identify individuals and matching the genotypes individuals' `IID`. In the examples, this file contains raw UK Biobank fields information (e.g. `1160-2.0`) that will be interpreted by the `UKB_CONFIG`.
- `UKB_WITHDRAWAL_LIST`: List of individuals to be removed from the analysis (one `IID` per line and no header).
- `BED_FILES`: This parameter points to the [plink bed files](https://www.cog-genomics.org/plink/1.9/) containing typed genetic variants information. You will notice the "{1,2,3}.{bed,bim,fam}" suffix of the parameter which is telling Nextflow to group chromosome files together.
- `BGEN_FILES`: Contain imputed variants information in [BGEN](https://www.chg.ox.ac.uk/~gav/bgen_format/) format. You will notice the "{1,2,3}.{bgen,bgen.bgi,sample}" suffix of the parameter which is telling Nextflow to group chromosome files together.

## 5. Analyse the Results

Hopefully, this is where you'll find something new or interesting in some way. This section is up to you but we recommend to read the section [Understanding TarGene's Outputs](@ref) or you will likely be a little lost...

!!! warning "Persistent Session"
    Most TarGene runs will not be rapid. This is because the pipeline includes multiple steps and some of them are computationally intensive. For that reason it is recommended to use a terminal multiplexer like [tmux](https://github.com/tmux/tmux/wiki) or [screen](https://www.gnu.org/software/screen/manual/screen.html).
