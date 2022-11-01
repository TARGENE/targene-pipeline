# Contributing

Contributions, whether bug fixes, new features or documentation improvements are very welcome!

## Raise an issue

In order to discuss and track the evolution of the project, please first raise an issue on the [targene-pipeline](https://github.com/TARGENE/targene-pipeline/issues) repository. If a change is agreed upon, the discussion should identify the relevant repositories that are concerned by the change. For instance, if one wishes to improve the extraction of traits from the UK-Biobank, the [UKBMain.jl](https://github.com/TARGENE/UKBMain.jl) would surely be impacted and a new release for that package necessary. For now, at the pipeline level, each release is considered "breaking". This means that to propagate the changes, a new release, updating the dependency on the UKBMain package must also be made for the targene-pipeline.

## Developping

For each repository identified in the previous step, create a new branch for your contribution. When you think your work is ready, open a pull request for review.

## Releasing

Apart from the targene-pipeline, a release consists in a docker image published on the docker hub. For each repository, by creating a github release with an appropriate tag, an action will be triggered and the image built and pushed to the registry.