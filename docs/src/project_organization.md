# Project Organization

The TarGene project is organized around the [targene-pipeline](https://github.com/TARGENE/targene-pipeline) repository which contains the [Nextflow](https://www.nextflow.io/) pipeline which will be the entry-point for most users. However, this repository does not contain the executables that are used by the Nextflow processes. Those executables originate from additional repositories:

- [UKBMain.jl](https://github.com/TARGENE/UKBMain.jl)
- [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl)
- [TargetedEstimation.jl](https://github.com/TARGENE/TargetedEstimation.jl)

The following diagram presents a high level perspective of the project's organization.

```@raw html
<div style="text-align:center">
<img src="assets/targene_organization.png" alt="TarGene Organization" style="width:800px;"/>
</div>
```

