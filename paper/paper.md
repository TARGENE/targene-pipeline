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
    affiliation: "1, 2"
  - name: Breeshey Roskams-Hieter
    affiliation: "1, 2"
    orcid: 0000-0002-1119-2576
  - name: Joshua Slaughter
    affiliation: "1, 2"
    orcid: 0000-0002-1400-6599
  - name: Kelsey Tetley-Campbell
    affiliation: "1, 2"
  - name: Mark J. van der Laan
    affiliation: "4"
  - name: Chris P. Ponting
    affiliation: "1"
    orcid: 0000-0003-0202-7816
  - name: Ava Khamseh
    affiliation: "1, 2, 4"
    orcid: 0000-0001-5203-2205
  - name: Sjoerd Viktor Beentjes
    affiliation: "1, 3, 4"
    orcid: 0000-0002-7998-4262
affiliations:
  - name: MRC Human Genetics Unit, Institute of Genetics and Cancer, University of Edinburgh, Edinburgh EH4 2XU, United Kingdom.
    index: 1
  - name: School of Informatics, University of Edinburgh, Edinburgh EH8 9AB, United Kingdom
    index: 2
  - name: School of Mathematics and Maxwell Institute for Mathematical Sciences, University of Edinburgh, Edinburgh EH9 3FD, United Kingdom
    index: 3
  - name: Division of Biostatistics, University of California, Berkeley, CA, USA
    index: 4

date: 14 January 2025
bibliography: paper.bib
---

# Summary

Genetic variants are the foundation of biological diversity, they play a crucial role in the adaptability, survival, and evolution of populations. Discovering which and how genetic variants affect human traits is an ongoing challenge with applications in healthcare and medicine. In some cases, genetic variants have an obvious effect because they change the coding sequence of a gene and thus its function. In the vast majority of cases however, variants occur in sequences of unknown function and could impact human traits or disease mechanisms in complex ways. TarGene is a Nextflow pipeline leveraging highly flexible machine-learning methods and semi-parametric estimation theory to capture these complex genetic dependencies including higher-order interactions.

# Statement of Need

All currently existing software for the estimation of genetic effects are based on parametric distributions, additionally assuming linearity of the relationship between variants and traits [@purcell2007plink,yang2011gcta,loh2018mixed,zhou2018efficiently]. If these assumptions are violated, the reported effect sizes will be biased and error rates inflated. In particular, this can lead to inflated false discovery rates and suboptimal allocation of computational resources and research funding. Some recently published software also account for more complex relationships but do not offer the full modelling flexibility provided by TarGene. REGENIE fits a two-stage whole-genome model for each phenotype of interest but still assumes linearity and normality [@mbatchou2021computationally]. DeepNull is a semi-parametric method which models non-linear covariate effects but also assumes genetic effects to be linear and does not allow complex interactions between covariates and genetic variants [@mccaw2022deepnull]. KnockoffGWAS [@sesia2021false] is non-parametric but does not estimate effect sizes, instead it aims at controlling the false discovery rate of variant selection in a genome-wide manner. In comparison, TarGene is the only method able to model arbitrarily complex genetic effects while preserving the validity of statistical inference. It does so by leveraging Targeted Learning [@van2011targeted], a framework combining methods from causal inference, machine learning and semi-parametric statistical theory. Succinctly, the estimation process works as follows. In a first step, flexible machine-learning algorithms are fitted to the data, hence minimizing an appropriate loss function (e.g., negative log-likelihood). A second step, known as the targeting step, regularises the estimate of the quantity of interest in a theoretically optimal way.

# Features

TarGene is a fully featured command-line software, which can be run as follow:

```
nextflow run https://github.com/TARGENE/targene-pipeline/ \
  -r TARGENE_VERSION \
  -c CONFIG_FILE \
  -resume
```

where the `CONFIG_FILE` provides the list of problem-specific parameters (data, arguments, options). Below we list some important features of TarGene, the following `CONFIG_FILE` will serve as a running example.

```
params {
    ESTIMANDS_CONFIG = "gwas_config.yaml"
    ESTIMATORS_CONFIG = "wtmle--tunedxgboost"

    // UK Biobank specific parameters
    BED_FILES = "unphased_bed/ukb_chr{1,2,3}.{bed,bim,fam}"
    UKB_CONFIG = "ukbconfig_gwas.yaml"
    TRAITS_DATASET = "dataset.csv"
}
```

For detailed explanations, please refer to the online [documentation](https://targene.github.io/targene-pipeline/stable/).

## Scalability

Machine learning methods are computationally intensive, however statistical genetics analyses need to scale to hundreds of thousands of variants and thousands of traits. For this reason, TarGene leverages Nextflow [@di2017nextflow], a pipeline management system that can parallelize independent estimation tasks across HPC platforms.

## Databases

TarGene works with standard formats, plink `.bed` and `.bgen` formats for genotypes, `.csv` or `.arrow` format for human traits. Furthermore, TarGene has direct support for two large scale biomedical databases, the UK Biobank [@bycroft2018uk] and the All of Us cohort [@all2019all]. The example considers the UK Biobank for which genotypes and traits are provided via `BED_FILES` and `TRAITS_DATASET` respectively. Because the UK Biobank has a non-standard format, the `UKB_CONFIG` provides trait definition rules. The following is an illustration for the body mass index phenotype, but the default is to consider all 768 traits as defined by geneAtlas [@canela2018atlas].

```
traits:
  - fields:
      - "21001"
    phenotypes:
      - name: "Body mass index (BMI)"
```

## Study Designs

TarGene supports traditional study designs in population genetics, that is, genome-wide association studies (GWAS) and phenome-wide association studies (PheWAS). Because TarGene has a focus on complex effects, interactions (e.g. gene-gene, gene-environment, gene-gene-environment) can also be investigated up to any order.

The study design is specified in the `ESTIMANDS_CONFIG` YAML file. For a routine GWAS the content of this file can be as simple as:

```
type: gwas
```

## Estimators

Semi-parametric estimators exist in multiple flavors, all with different properties. In TarGene we default to using Targeted Maximum-Likelihood Estimation [@van2018targeted] and XGBoost [@chen2016xgboost] as the machine-learning model. We have selected this default because it was the best performing estimator in simulations across a variety of genetics tasks [@labayle2025semi]. In the presence of computational restrictions, tradeoffs can be made and lighter models can be used.

# Acknowledgements

Olivier Labayle was supported by the United Kingdom Research and Innovation (grant EP/S02431X/1), UKRI Centre for Doctoral Training in Biomedical AI at the University of Edinburgh, School of Informatics.
Breeshey Roskams-Hieter was supported by the Health Data Research UK & The Alan Turing Institute Wellcome PhD Programme in Health Data Science (Grant Ref: 218529/Z/19/Z).
Mark van der Laan is supported by NIH grant R01AI074345.
Chris P. Ponting was funded by the MRC (MC_UU_00007/15).
Ava Khamseh was supported by the XDF Programme from the University of Edinburgh and Medical Research Council (MC_UU_00009/2), and by a Langmuir Talent Development Fellowship from the Institute of Genetics and Cancer, and a philanthropic donation from Hugh and Josseline Langmuir.

# References {#references .unnumbered}
