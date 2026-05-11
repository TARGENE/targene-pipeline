# TarGene-Pipeline

TarGene is a Nextflow pipeline for the estimation of genetic effects on human traits using highly flexible semi-parametric inference methods (Targeted Minimum Loss-Based Estimation, One-Step Estimation). 

## Is TarGene for you?

If you would like to use flexible machine-learning methods while preserving valid statistical inference **and** your project falls within any of the following, then the answer is Yes!

Supported study designs:

- GWAS
- PheWAS
- Custom

Supported databases:

- UK Biobank
- All of Us
- Custom

Supported genetic effects:
- Main effects
- Interactions up to any order

## Documentation Access

Want to know more? Please visit [the docs](https://targene.github.io/targene-pipeline/stable/).

## Running the Workflow

TarGene is a Nextflow pipeline, please refer to [their documentation](https://www.nextflow.io/) for general usage. The main point being that, depending on your HPC specifications, you will need to provide a specific `myplatform.config` configuration file on top of your project `myproject.config` configuration file.

```bash
nextflow run TARGENE/targene-pipeline -c myplatform.config -c myproject.config -resume
```
