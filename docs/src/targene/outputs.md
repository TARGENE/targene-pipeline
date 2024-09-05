# Understanding TarGene's Outputs

The successful completion of the workflow will produce the following files in the output directory `OUTDIR (default: results)`:

- `QQ.png`: A Quantile-Quantile summary plot.
- `results.summary.yaml`: A summary statistics file of your results.
- `results.hdf5`: A file containing the complete description of estimands and associated estimates in HDF5 format.
- `svp.hdf5`: An optional file containing Sieve Variance Plateau corrected variance estimates (see [Correcting for population relatedness](@ref)).

If you come from the world of linear models where a single ``\beta`` and p-value is output for each variant, TarGene's output may seem difficult to read at first. In this section we explain how they are structured and how to work with them. 

The crucial difference, is that TarGene estimates the effect of change. As such there is one estimate (``\hat{\beta}``) per genotype change. Consider a single variant with genotypes CC, CT and TT. TarGene estimates non-redundant changes, in this case it would be both ``CC \rightarrow CT`` and ``CT \rightarrow TT``. For each of these changes, there is thus an effect size, e.g. ``\hat{\beta}_{CC \rightarrow CT}``, and a p-value, e.g. ``\text{pval}_{\hat{\beta}_{CC \rightarrow CT}}``. 

However, it is often more intuitive to think of the "full" effect of a variant on an outcome. This is simply defined as the joint effect of all changes, i.e., ``[\beta_{CC \rightarrow CT}, \beta_{CT \rightarrow TT}]``. This effect is nonzero if it is different from ``[0, 0]``. The multivariate central limit theorem tells us that this joint effect is a multivariate normal and a single p-value can be obtained for it.

Finally, estimation can only be accurate for genotypes that are not too rare. This is known as the [positivity condition](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8492528/), and is similar to the traditional minor allele frequency threshold used in GWAS. In TarGene this is defined by the `POSITIVITY_CONSTRAINT (default: 0.01)` parameter, for which a default value was chosen based on simulation studies. In TarGene, only changes that satisfy the positivity constraint are considered for estimation. In the example above, if the genotypes frequencies were 0.745, 0.25 and 0.005 respectively, only the ``CC \rightarrow CT`` would satisfy the constraint, and only this change would be estimated.

We have discussed single variant effects, but interactions are defined exactly in the same way. The difference is that they are defined by two or more variants, each of which with possibly multiple genotype changes. Consider two variants ``V_1`` and ``V_2``, the interaction between ``V_1`` and ``V_2`` is defined by the following genotype changes:

- ``V_1: CC \rightarrow CT, V_2: AA \rightarrow AG``
- ``V_1: CT \rightarrow TT, V_2: AA \rightarrow AG``
- ``V_1: CC \rightarrow CT, V_2: AG \rightarrow GG``
- ``V_1: CT \rightarrow TT, V_2: AG \rightarrow GG``


## The `results.summary.yaml`

## The `results.hdf5`

