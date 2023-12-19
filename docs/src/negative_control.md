# Running negative control checks

After running TarGene, you may want to perform some quality control checks on your positive hits. Since true positives are not ubiquitous in biology, this is done by negative control checks. The corresponding pipeline can be run by:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -entry negativeControl -r vX -profile P -resume -with-trace -with-report
```

see the [Overview](@ref) section for more info on the command line options.

At the moment, there are two types of analyses you can perform with this sub-workflow: permutation tests and comparisons with "non-functional" random variants drawn from the genome. In both cases, you will need to define positive hits. This is traditionally done via p-value filtering. Since TarGene allows the use of multiple estimators, you need to specify for which one you would like to run negative controls. This is done via the `ESTIMATOR_KEY` (default: "TMLE") parameter. You can also modify the `PVAL_THRESHOLD` (default: 0.05), only those estimates below the threshold are considered for negative control.

## Permutation tests

For each IATE hit, we can perturb the data by independently shuffling one or more columns in the data. By doing so, the interaction will be broken and become "non-significant". Because there are up to K variables in a IATE estimand, we can also permute up to K+1 variables. For instance, assuming we found a significant interaction between rs1234 and rs5678 on diabetes, we could permute:

- Any one of (rs1234, rs5678, diabetes) leading to 3 order 1 permutation tests.
- Any two of (rs1234, rs5678, diabetes) resulting also in 3 order 2 permutation tests.
- All three of (rs1234, rs5678, diabetes) leading to 1 order 3 permutation test.

Specifying those orders is made via the `PERMUTATION_ORDERS` variable which is a comma separated string of orders (e.g. "1,2").

Just like the main pipeline, the outputs will be found in:

- `PERMUTATION_HDF5_OUTPUT` (default: "permutation_results.hdf5")
- `PERMUTATION_JSON_OUTPUT` if specified

## Non-functional randomly drawn variants

If you are using the `STUDY_DESIGN` = `FROM_ACTORS`, it is likely that you have not chosen the genetic variants at random. Mainly trans-acting variants have been defined because of a likely role in a biological mechanism. If this is the case, replacing a trans-actor by a random variant anywhere on the genome is likely to break the interaction. While unlikely, by chance alone (or lack thereof), one could pick a random trans-actor which is also interacting with any of the bQTLs. In order to account for that it is recommended to instead select a certain number of random variants, denoted by `N_RANDOM_VARIANTS` (default: 10). Furthermore, we enforce two criteria on each of the randomly chosen random variants:

- Its minor allele frequency should match that of the original trans-actor up to `MAF_MATCHING_RELTOL` (relative tolerance, default: 0.05).
- It shouldn't lie in a known regulatory region.

The result of this part of the pipeline is a file named `random_variants_estimands.jls`. To perform the actual tests, you will have to run TarGene again on that using this file with the `STUDY_DESIGN` = `CUSTOM` mode.
