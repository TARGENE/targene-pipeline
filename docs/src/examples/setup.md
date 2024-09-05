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

Make sure the data is available. Here I will rely on TarGene's test data located [here](https://github.com/TARGENE/targene-pipeline/tree/main/test/assets). To obtain it, the easiest way is maybe to clone the entire project like so 

```bash
git clone https://github.com/TARGENE/targene-pipeline
```

You can then open the `test/assets` directory (or simply the full project but you need to open a terminal in `test/assets`) using your favorite text editor. 

!!! note "This is only to get the data"
    Even though we have cloned the repository, we do not need to look at the code at all!

All examples below assume that the dataset is the UK-Biobank data (we are not really using the UK-Biobank but the data was created to mimic its structure).

## 3. Writing a Run Configuration File

Write a run configuration file for your project, this can be decomposed further into 2 sub-steps:

1. Write the part of your configuration file specific to your platform. Here we will simply use the `local` profile.
2. Write the part of your configuration file specific to your run. This is the topic of the remainder of the following sections.

## 4. Analyse the Results

Hopefully, this is where you'll find something new or interesting in some way. This section is up to you but we recommend to read the section [Understanding TarGene's Outputs](@ref) or you will likely be a little lost...
