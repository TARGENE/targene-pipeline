---
title: 'TarGene: A Nextflow pipeline for the estimation of genetic effects on human traits via semi-parametric methods.'
tags:
  - nextflow
  - statistics
  - semi-parametric statistics
  - population genetics
  - genetic effects
  - causal inference
authors:
  - name: Olivier Labayle
    orcid: 0000-0002-3708-3706
    equal-contrib: true
    affiliation: "1, 2"
  - name: Breeshey Roskams-Hieter
    affiliation: "1, 2"
    orcid: 0000-0002-1119-2576
  - name: Joshua Slaughter
    affiliation: "1, 2"
    orcid: 0000-0002-1400-6599
  - name: Kelsey Tetley-Campbell
    affiliation: "1, 2"
  - name: Chris P. Ponting
    affiliation: "3"
    orcid: 0000-0003-0202-7816
  - name: Sjoerd Viktor Beentjes
    affiliation: "1, 3, 4"
    orcid: 0000-0002-7998-4262
  - name: Ava Khamseh
    affiliation: "1, 2, 4"
    orcid: 0000-0001-5203-2205
affiliations:
  - name: MRC Human Genetics Unit, Institute of Genetics and Cancer, University of Edinburgh, Edinburgh EH4 2XU, United Kingdom.
    index: 1
  - name: School of Informatics, University of Edinburgh, Edinburgh EH8 9AB, United Kingdom
    index: 2
  - name: School of Mathematics and Maxwell Institute, University of Edinburgh, Edinburgh EH9 3FD, United Kingdom
    index: 3
  - name: Division of Biostatistics, University of California, Berkeley, CA, USA
    index: 4

date: 14 January 2025
bibliography: paper.bib
---

# Summary

Genetic variations are the foundation of biological diversity, they play a crucial role in the adaptability, survival, and evolution of populations. Discovering which and how genetic variations affect human traits is an ongoing challenge with applications in healthcare and medicine. In some cases, genetic variations have an obvious effect because they change the coding sequence of a gene and thus its function. In the vast majority of cases however, variations occur in places of unknown function and could impact human traits or disease mechanisms in complex ways. TarGene is a Nextflow pipeline leveraging highly flexible machine-learning methods and semi-parametric estimation theory to capture these complex genetic dependencies including higher-order interactions.

# Statement of Need

All currently existing software for the estimation of genetic effects are based on parametric models, usually assuming linearity and normality. If these assumptions are violated, the effect sizes reported by these sofware will be biased. This can lead to spurious associations, hence inflating false discovery rates and leading to suboptimal resources allocation. We list here a few recent methods that address some of the limitations pointed out. REGENIE has the benefit to fit a whole-genome model for each phenotype of interest but still assumes linearity and normality [@mbatchou2021computationally]. DeepNull is a semi-parametric method which models non-linear covariate effects but still assumes genetic effects to be linear and does not allow complex interactions between covariates and genetic variants [@mccaw2022deepnull]. KnockoffGWAS, aims at controlling the false discovery rate in genome-wide association studies. It does not rely on strong parametric assumptions but does not estimate effect sizes [@sesia2021false].

# Features

## Nextflow

Machine-learning methods are computationally more intensive than traditional linear models. For this reason, TarGene leverages Nextflow [@di2017nextflow], a pipeline management system that can parallelize independent estimation tasks across HPC platforms.

TarGene can be run via the command line

```
nextflow run https://github.com/TARGENE/targene-pipeline/ -r TARGENE_VERSION -c CONFIG_FILE -resume
```

Below we list some important features of TarGene. The following `CONFIG_FILE` will serve as a running example.

```
params {
    ESTIMANDS_CONFIG = "gwas_config.yaml"
    ESTIMATORS_CONFIG = "wtmle--glm"

    // UK-Biobank specific parameters
    BED_FILES = "unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    UKB_CONFIG = "ukbconfig_gwas.yaml"
    TRAITS_DATASET = "dataset.csv"
}
```

## Databases

TarGene works with standard formats, plink `.bed` and `.bgen` formats for genotypes, `.csv` or `.arrow` format for human traits. Furthermore, TarGene has direct support for two large scale biomedical databases, the UK-Biobank and the All of Us cohort. The example considers the UK-Biobank for which genotypes and traits are provided via `BED_FILES` and `TRAITS_DATASET` respectively. Because the UK-Biobank has a non-standard format, the `UKB_CONFIG` provides traits definition rules. The following is an illustration for BMI, but the default is to consider all 766 traits as defined by the geneAtlas [@canela2018atlas].

```
traits:
  - fields:
      - "21001"
    phenotypes:
      - name: "Body mass index (BMI)"
```

## Study Designs

TarGene supports traditional study designs in population genetics, that is, genome-wide association studies (GWAS) and phenome-wide association studies (PheWAS). Because TarGene has a focus on complex effects, higher-order interactions (e.g. gene-gene-... or gene-environment-...) can also be investigated.

The study design is specified in the `ESTIMANDS_CONFIG` YAML file. For a routine GWAS this can be as simple as:

```
type: gwas
```

## Estimators

Semi-parametric estimators exist in multiple flavors, all with different properties. In TarGene we default to using Targeted Maximum-Likelihood Estimation [@van2018targeted] and XGboost [@chen2016xgboost] as the machine-learning model. This is because this was the best performing estimator in simulations for a variety of tasks. But if computational restrictions exist, tradeoffs must be made. In the running example, XGBoost is replaced by a generalised linear model.

# Mention

TarGene
- was used to evaluate the performance of semi-parametric estimators in real-world data scenarios using the the UK-Biobank in
- is being used to discover causal genetic variants acting through biological mechanisms via interaction analyses 


# Acknowledgements

This work was supported by the United Kingdom Research and Innovation (grant EP/S02431X/1), UKRI Centre for Doctoral Training in Biomedical AI at the University of Edinburgh, School of Informatics.

# References {#references .unnumbered}
