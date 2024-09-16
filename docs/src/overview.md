# Nextflow Basics

## Running the Workflows

Since TarGene uses [Nextflow](https://www.nextflow.io/), all workflows can be run in the same way from the command line:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline/ -r TARGENE_VERSION -entry WORKFLOW_NAME -profile P -resume
```

where:

- `TARGENE_VERSION` is the latest TarGene version, e.g. `v0.11.0`
- `WORKFLOW_NAME` is any of the [TarGene workflows](@ref "Project Configuration")
- `P` is an optional [Nextflow profile](https://www.nextflow.io/docs/latest/config.html) describing the computing platform (see [Platform Configuration](@ref)).

Additional Nextflow command line arguments can be found in the official documentation, for example important options are:

- `-resume`: Tells Nextflow to try to resume the pipeline if an error occurred during the execution (if you forgot to specify a parameter for instance)
- `-with-trace` and `-with-report` will generate additional report files.

## Workflows Configurations

There are mainly two parts to configuring a workflow run. The first part describes the computing environment, it tells Nextflow how it should execute the various processes of the workflow. The second part provides the actual inputs to the TarGene workflows. When running the `nextflow run` command, Nextflow will look for a `nextflow.config` configuration file in your current directory. If you are new to Nextflow, you can use this file to setup both the platform and project configurations. As your project grows you may want to split them in [distinct files](https://www.nextflow.io/docs/latest/config.html#configuration). This is so that the platform configuration can be easily reused across many TarGene runs.

### Platform Configuration

It is likely that you will run TarGene on a HPC platform, in particular the [Executors](https://www.nextflow.io/docs/latest/executor.html) and [Singularity](https://www.nextflow.io/docs/latest/container.html#singularity) configurations are required. Since Nextflow is so widespread, it is probable that such a configuration file is already available from your HPC administrators. Since this configuration only describes de computing platform and not your project, it is often described as a [Profile](https://www.nextflow.io/docs/latest/config.html#config-profiles).

!!! tip "University of Edinburgh"
    If you are using the University of Edinburgh [Eddie cluster](https://www.ed.ac.uk/information-services/research-support/research-computing/ecdf/high-performance-computing), you can simply use TarGene with the `-profile eddie` option.

### Project Configuration

These are the configuration details associated with your project, this is usually done in a `nextflow.config` file living at the root of your project's directory. The configuration parameters are specific to each workflow and described in the following sections.
