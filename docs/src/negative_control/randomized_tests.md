# The Randomized Variants Workflow

If you are using the `STUDY_DESIGN` = `FROM_ACTORS`, it is likely that you have not chosen the genetic variants at random. Mainly trans-acting variants have been defined because of a likely role in a biological mechanism. If this is the case, replacing a trans-actor by a random variant anywhere on the genome is likely to break the interaction. While unlikely, by chance alone (or lack thereof), one could pick a random trans-actor which is also interacting with any of the bQTLs. In order to account for that it is recommended to instead select a certain number of random variants, denoted by `N_RANDOM_VARIANTS` (default: 10). Furthermore, we enforce two criteria on each of the randomly chosen random variants:

- Its minor allele frequency should match that of the original trans-actor up to `MAF_MATCHING_RELTOL` (relative tolerance, default: 0.05).
- It shouldn't lie in a known regulatory region.

The result of this part of the pipeline is a file named `random_variants_estimands.jls`. To perform the actual tests, you will have to run TarGene again on that using this file with the `STUDY_DESIGN` = `CUSTOM` mode.

## Example Run Command

```bash
nextflow run https://github.com/TARGENE/targene-pipeline/ -r TAG -entry RANDOMIZATION_TEST -profile P -resume
```

## List Of Workflow Arguments

- **`NB_PCS` (optional, default: 6)**: The number of PCA components to extract.
- **`BED_FILES` (required)**: Path expression to PLINK BED files.
- **`COHORT` (optional: "UKBB")**: Current default for this is UKBB. If set to a value other than UKBB, this will not run UKBB-specific trait extraction.
- **`TRAITS_DATASET` (required)**: Path to a traits dataset. If you are running this for a non-UKBB cohort, your sample IDs must be specified in the first column of this CSV file, with the column name `SAMPLE_ID`.
- **`FLASHPCA_EXCLUSION_REGIONS` (optional, default: data/exclusion_regions_hg19.txt)**: A path to the flashpca special exclusion regions.
- **`MAF_THRESHOLD` (optional, default: 0.01)**: Only variants with that minor allele frequency are considered
- **`LD_BLOCKS` (optional)**: A path to pre-identified linkage disequlibrium blocks to be removed from the BED files. It is good practice to specify `LD_BLOCKS`, as it will remove SNPs correlated with your variants-of-interest before running PCA.

**If the `COHORT` argument is set to `UKBB`**:

- **`UKB_CONFIG` (required)**: YAML configuration file describing which traits should be extracted and how the population should be subsetted.
- **`UKB_ENCODING_FILE` (optional)**: If the `TRAITS_DATASET` is encrypted, an encoding file must be provided.
- **`UKB_WITHDRAWAL_LIST` (optional)**: List of participants withdrawn from the study.
- **`QC_FILE` (optional)**: Genotyping quality control file from the UK-Biobank study.
