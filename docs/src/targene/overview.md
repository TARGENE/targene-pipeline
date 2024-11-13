# Overview

## General Workflow Structure

This is the main workflow within TarGene, its purpose is to estimate a wide variety of genetic effects using the Targeted Learning framework. This is an end-to-end workflow, meaning that you don't need to perform any QC on your genotypes files. The workflow can be roughly decomposed into two main steps:

1. In the first step, an integrated tabular dataset is built, it contains
   - Phenotypes: Potentially extracted from the UK Biobank
   - Variants of Interest: Extracted from genotyping data.
   - PCs: Constructed from genotyping data using standard methodology (A LOCO approach is used for GWAS)
2. In the second step, all genetic effects are estimated via Targeted Learning in parallel using the estimators of your choice.

An overview of the workflow is presented in the following diagram.

![TarGene Workflow High Level](../assets/targene_workflow_high_level.png)

## Example Run Command

```bash
nextflow run https://github.com/TARGENE/targene-pipeline/ -r v0.11.1 -profile local -resume
```

We now describe step by step how to setup a TarGene run configuration.
