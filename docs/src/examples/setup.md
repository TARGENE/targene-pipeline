# Setup: Read First

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

Make sure the data is available. Here I will rely on TarGene's test data located [here](https://github.com/TARGENE/targene-pipeline/tree/main/test/assets). To obtain it, the easiest way is maybe to clone the project like so 

```bash
git clone https://github.com/TARGENE/targene-pipeline
```

and open the `test/assets` directory.

All examples below assume that the dataset is the UK-Biobank data (we are not really using the UK-Biobank but the data was created to mimic its structure). We will use:

- `TRAITS_DATASET` (e.g. `dataset.csv`): A raw dataset as downloaded from the UK-Biobank's portal and containing fields data. Here it is assumed to have been decrypted.
- `UKB_CONFIG` (e.g. `ukbconfig_gwas.yaml`): A list of traits and phenotypes to be extracted from the `TRAITS_DATASET`.
- `UKB_WITHDRAWAL_LIST` (e.g. `withdrawal_list.txt`): A list of participants that have been withdrawn from the study.
- `BED_FILES` (e.g. `unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}`): A list of unphased [plink BED files](https://zzz.bwh.harvard.edu/plink/binary.shtml) containing the genotyping data.

## 3. Writing a Run Configuration File

Write a run configuration file for your project, this can be decomposed further into 2 sub-steps:

1. Write the part of your configuration file specific to your platform. Here we will simply use the `local` profile.
2. Write the part of your configuration file specific to your run. This is the topic of the remainder of the following sections.

## 4. Analyse the Results

Hopefully, this is where you'll find something new or interesting in some way. This section is up to you but we recommend to read the section [Understanding TarGene's Outputs](@ref) or you will likely be a little lost...
