# Understanding TarGene's Outputs

The successful completion of the workflow will produce the following files in the output directory `OUTDIR (default: results)`:

- `QQ.png`: A Quantile-Quantile summary plot:
- `results.summary.yaml`: A summary file of your results.
- `results.hdf5`: A file containing the complete description of estimands and estimators in HDF5 format.
- `svp.hdf5`: An optional file containing Sieve Variance Plateau corrected variance estimates (see [Correcting for population relatedness](@ref)).

If you come from the world of linear models where a simple ``\beta`` and p-value is output for each variant, TarGene's output may seem difficult to read at first. In this section we explain how they are structured and how to work with them. 

The crucial difference is that TarGene estimates the effect of change. As such there is one estimate (``\beta``) per genotype change. Consider a single variant with genotypes CC, CT and TT. TarGene estimates non-redundant changes, in this case it would be both ``CC \rightarrow CT`` and ``CT \rightarrow TT``. For each change there is thus an effect size, e.g. ``\hat{\beta}_{CC \rightarrow CT}``, and a p-value, e.g. ``\text{pval}_{\hat{\beta}_{CC \rightarrow CT}}``. 

However, it is often more intuitive to think of the "full" effect of a variant on an outcome. We say that this is the case if the joint effect ``[\beta_{CC \rightarrow CT}, \beta_{CT \rightarrow TT}]`` is different from ``[0, 0]``. Since this joint effect is a multivariate normal, this can be tested and a single p-value can be obtained for it.

## The `results.summary.yaml`

## The `results.hdf5`

