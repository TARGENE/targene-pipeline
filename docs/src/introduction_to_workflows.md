# Workflows

There are currently 2 main workflows and two secondary workflows within TarGene. Since TarGene uses [Nextflow](https://www.nextflow.io/), all workflows can be run in the same way from the command line:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline/ -r TAG -entry WORKFLOW_NAME
```

where `TAG` is the latest TarGene version, e.g. `v0.9.0` and `WORKFLOW_NAME` is any of the worflows listed below.

## Main Workflows

- [The TarGene Workflow](@ref) (`WORKFLOW_NAME: TARGENE`): It is the main workflow for the targeted estimation of genetic effects.
- [The Negative Control Workflows](@ref): These workflow enable the control of the false discovery rate by using the results obtained from a previous TarGene discovery run. There are currently two of them:
   - [Permutation Test Workflow](@ref) (`WORKFLOW_NAME: PERMUTATION_TEST`): Performs permutation tests by independently shuffling the individuals in the columns of an aggregated dataset.
   - [Randomized Variants Workflow](@ref) (`WORKFLOW_NAME: RANDOMIZATION_TEST`): When there are multiple genetic variants of interest, e.g. in an interaction study, one can replace one of the variant at random by another variant and the effect is expected to be 0 in average.

## Secondary Workflows

- [The PCA Workflow](@ref) (`WORKFLOW_NAME: PCA`): This workflow computes principal components.
- [The Make Dataset Workflow](@ref) (`WORKFLOW_NAME: MAKE_DATASET`): This workflow generates an aggregated dataset from traits and genetic data.

