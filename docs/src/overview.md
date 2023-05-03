# Overview

Since TarGene is a Nextflow pipeline, it can be run with a simple command line similar to the following:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r vX -profile P -resume -with-trace -with-report
```

All arguments are optional but encouraged. Here `-r vX` describes the version to be used and should be provided for reproducibility purposes, e.g. `-r v0.3.7`. If left unspecified, the latest development version will be used. The `-resume` option tells Nextflow to try to resume the pipeline if an error occured during the execution (if you misspecified a parameter for instance) and the `-with-trace` and `-with-report` generate additional report files. The `-profile P` option is described below and implicit, is the existence of a `nextflow.config` file, the content of which is also described here:

1. It is likely that you will run TarGene on a HPC platform, in particular the [Executors](https://www.nextflow.io/docs/latest/executor.html) and [Singularity](https://www.nextflow.io/docs/latest/container.html#singularity) configurations are required. Since Nextflow is so widespread, it is probable that such a configuration file is already available from your HPC administrators. Since this configuration only describes de computing platform and not your project, it is often described as a [Profile](https://www.nextflow.io/docs/latest/config.html#config-profiles). If your HPC uses the SGE executor, the `-profile eddie` may work with no, or minor adjustment (it can also serve as a template for other executors [see file](https://github.com/TARGENE/targene-pipeline/blob/main/conf/eddie.config)).

2. You need to provide the configuration details associated with your project, this is usually done in a `nextflow.config` file living at the root of your project's directory. The configuration parameters are described in the following sections:
    - [Setting a data source](@ref)
    - [Adjusting for confounders](@ref)
    - [Describing the causal parameters of interest](@ref)
    - [Specifying a Targeted Estimator](@ref)
    - [Correcting for population relatedness](@ref)
    - [Tweaking additional behaviour](@ref)
    - [Running negative control checks](@ref)

Finally, a list of all TarGene's parameters is available in the [Index of the pipeline parameters](@ref).
