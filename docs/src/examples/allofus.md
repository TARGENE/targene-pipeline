# All of Us

Analyses within the All of Us (AoU) Researcher Workbench using genetic data must be run within the `Controlled Tier Access` (See [`Data Access`](https://www.researchallofus.org/data-tools/data-access/)). Workspaces launched within this tier will automatically have nextflow installed and can use TarGene immediately. Each project will be assigned a bucket for storage on Google Cloud (this can be viewed ), as well as all relevant nextflow parameters required to submit jobs using the Google Lifesciences API. Your workspace-specific nextflow profile will be found in `~/.nextflow/config` and will be automatically available to you by using the flag `-profile gls` when you run nextflow. The AoU Researcher workbench requires some additional specifications that are built into TarGene in the `allofus` profile, and can be combined with your Workspace-specific `gls` configuration when running TarGene.

We reccommend running this by first entering a `Cloud Analysis Terminal` on your current Workspace, creating a configuration for the analysis you would like to run, and running TarGene in a screen session. See [`Workflows in the All of Us Researched Workbench`](https://support.researchallofus.org/hc/en-us/articles/4811899197076-Workflows-in-the-All-of-Us-Researcher-Workbench-Nextflow-and-Cromwell) for more information.

A minimalist run configuration to run a flat config run on the AoU Researcher Workbench using TarGene looks like the following:

```conf
params {
    COHORT = "ALLOFUS"
    ESTIMANDS_CONFIG = "allofus_config.yaml"

    // UK-Biobank specific parameters
    BGEN_FILE = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/clinvar_v7.1/bgen/clinvar.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}"
    BED_FILES = "gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/clinvar_v7.1/plink_bed/clinvar.chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}"
    TRAITS_DATASET = "allofus_traits.csv"
}
```

Apart from the data related parameters, there are two main parameters here: `ESTIMANDS_CONFIG` and the `ESTIMATORS_CONFIG`. These parameters describe the estimands (questions of interest) and how to estimate them respectively. 

The `ESTIMANDS_CONFIG` here follows the same format as the `flat` configuration detailed in the `PheWas` section. Here we are estimating the ATE of the FTO variant (encoded as chr16:53767042:T:C in the BGEN data for the AoU cohort) on the traits present in our `allofus_traits.csv` file. We have also added `sex_at_birth` as a covariate. 

```yaml
type: flat

estimands:
  - type: ATE

variants:
  - chr16:53767042:T:C

outcome_extra_covariates:
  - "sex_at_birth"
```

The optional `outcome_extra_covariates` are variables to be used as extra predictors of the outcome (but not as confounders). This information must be contained in the `allofus_traits.csv` file, along with your outcomes-of-interest. 

The `allofus_traits.csv` might look as follows:

| SAMPLE_ID | sex_at_birth | Height (cm) |
|-----------|--------------|-------------|
| 100000    | Male         | 180         |
| 100002    | Female       | 165         |
| 100004    | Male         | 175         |
| 100010    | Female       | 160         |
| ...       | ...          | ...         |

Here we have not specified any value for `ESTIMATORS_CONFIG`, and so the default, `ESTIMATORS_CONFIG = "wtmle-ose--tunedxgboost"`, will be used. This defines the estimation strategy, and more specifically, that we will be using Targeted Minimum-Loss Estimator as well as a One Step Estimator with a tuned XGBoost model to learn the outcome models (``Q_Y``) and propensity scores (``G``).

The TarGene can then be run on the AoU Researcher Workbench as follows:

```bash
nextflow run https://github.com/TARGENE/targene-pipeline -r v0.11.1 -profile gls,allofus
```