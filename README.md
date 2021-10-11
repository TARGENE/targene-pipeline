# UKBBEpistasisPipeline

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://olivierlabayle.github.io/UKBBEpistasisPipeline.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://olivierlabayle.github.io/UKBBEpistasisPipeline.jl/dev)
[![Build Status](https://github.com/olivierlabayle/UKBBEpistasisPipeline.jl/workflows/CI/badge.svg)](https://github.com/olivierlabayle/UKBBEpistasisPipeline.jl/actions)
[![Coverage](https://codecov.io/gh/olivierlabayle/UKBBEpistasisPipeline.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/olivierlabayle/UKBBEpistasisPipeline.jl)

Here we provide the main workflow for estimating genetic variants interactions responsible of traits in the UK Biobank. For that purpose, we rely on [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html) and [NextFlow](https://www.nextflow.io/).

## Running the pipeline on Eddie

According to the [Eddie specific Nextflow documentation](https://www.wiki.ed.ac.uk/display/ResearchServices/Bioinformatics), to run the workflow with the basic configuration using Singularity, from the main directory:

```bash
nextflow run -c /exports/igmm/eddie/BioinformaticsResources/nextflow/config/eddie.config main.nf
```

or via gitlab (if your $HOME/.nextflow/scm is correctly set):

```bash
nextflow run -c /exports/igmm/eddie/BioinformaticsResources/nextflow/config/eddie.config tfomics/EstimationPipeline -hub uoegitlab
```

To adapt the configuration to your specific platform or need, you will need to adapt the configuration.