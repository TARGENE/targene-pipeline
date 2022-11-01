# TL-Pipeline

> **Warning**:
>   This is rapidly evolving software that should be considered unstable.

Here we provide the main workflow for estimating genetic variants interactions responsible of traits in the UK Biobank using the Targeted Learning framework. For that purpose, we rely on [NextFlow](https://www.nextflow.io/), a software that helps the development of complex parallel and reactive workflows on clouds and clusters.

## Running the Workflow

Please refer to the main [NextFlow](https://www.nextflow.io/) documentation for general usage. The main point being that, depending on your cluster specifications, you will need to provide a specific `myprofile` configuration file. If you are part of the University of Edinburgh and simply using Eddie, then the `eddie` profile is already defined. Then simply run:

```bash
nextflow run TARGENE/targene-pipeline -profile myprofile -resume
```
