# Contributing

Contributions, whether bug fixes, new features or documentation improvements are very welcome!

## Raise an issue

In order to discuss and track the evolution of the project, please first raise an issue on the [targene-pipeline](https://github.com/TARGENE/targene-pipeline/issues) repository. If a change is agreed upon, the discussion should identify the relevant repositories that are concerned by the change and open an issue on each of the repository. For instance, if one wishes to improve the extraction of traits from the UK-Biobank, the [UKBMain.jl](https://github.com/TARGENE/UKBMain.jl) would surely be impacted and a new release for that package necessary.

## Workflow

Following our previous UKBMain.jl example, there are two repositories tha need to be updated, the current workflow is as follows:

1. Develop
    - UKBMain.jl
        - Create a new git branch for your change
        - Develop and test
        - Release an image for your branch by selecting it after clicking the [Run workflow button](https://github.com/TARGENE/UKBMain.jl/actions/workflows/Release.yml). If the tests pass, a new docker image will be generated and hosted on Dockerhub with your branch's name.
    - targene-pipeline
        - Create a new git branch for your change
        - For each Nextflow process using the UKBMain.jl's docker image, update to the branch's image name.
        - Develop further required changes and run/add the tests (see [Note on the pipeline's tests](@ref)).
2. Review: When everything is working, ask for a review
3. Release UKBMain.jl:
    - Merge your branch into main
    - Create a new Github release, this will create a new docker image with your release name.
4. Release TarGene
    - For each Nextflow process using the UKBMain.jl's docker image, update to the released image name (as before).
    - Create a new Github release

## Note on the pipeline's tests

Currently the pipeline's end to end tests are not automated, which means you have to run them manually. Furthermore they only check that the pipeline terminates. We have a suite of two tests that can be run on Eddie by the following:

```bash
nextflow run main.nf -c conf/ci_jobs/from_actors.config -profile eddie -resume
```

and

```bash
nextflow run main.nf -c conf/ci_jobs/from_param_files.config -profile eddie -resume
```

Further functionalities should be accompanied with further test configurations.
