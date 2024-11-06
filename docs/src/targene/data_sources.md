# Genetic and Traits Data

The dataset is the main ingredient manipulated by TarGene. Population Genetics datasets are typically split into a variety of files and formats. The following section describes these formats and how to input them to a TarGene run.

## Genetic Data

TarGene currently only works with genotyping arrays datasets. These datasets can come in multiple formats and we currently use both `bgen` and `bed` formats. It is also assumed that the genotyping data is **unphased**. and split in chromosomes.

Unless you are performing a GWAS, both `bgen` and `bed` files are required. This is because PCA uses the `bed` files while the variants of interest are looked up within the more comprehensive `bgen` files. Since the later format does not provide variant calls but only imputation probabilities, the `CALL_THRESHOLD (default: 0.9)` parameter is used to call genotypes. In the GWAS setting, the complete list of all variants within the `.bed` files is used for analysis.

The `bgen` and `bed` files are provided via the `BGEN_FILES` and `BED_FILES` Nextflow parameters respectively. For example:

- `BGEN_FILES="path_to_folder/imputed_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}"` 
- `BED_FILES="path_to_folder/genotyped_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}`

!!! tip "File Names"
    It is advised to deviate as little as possible from the prefix above (prefix ending in "chr") to make sure everything works smoothly.

## Traits Dataset

In addition to genetic data, a set of traits is usually also available in a different format. We currently support two types of cohorts: the UK Biobank and custom CSV datasets.

### Custom Traits Dataset

Select with `COHORT = "MY_COHORT"` (The actual name does not matter as long as it does not match an existing cohort).

TarGene supports custom traits datasets in CSV format, that can be provided via the `TRAITS_DATASET` parameter. Please ensure the sample-ids identifying individuals in your cohort are included as the first column of your trait data, with the column name `SAMPLE_ID`.

### UK-Biobank

Select with `COHORT = "UKB"`.

If the dataset provided as `TRAITS_DATASET` is the encrypted dataset (potentially suffixed by `enc_ukb`), then the decryption file must be provided via `UKB_ENCODING_FILE`. If you have already decrypted the dataset, you don't need to provide a `UKB_ENCODING_FILE`. We refer to either the encrypted or decrypted dataset as the main dataset.

The information contained in the UK-Biobank main dataset is not directly usable for statistical analyses and needs to be processed. This is both because the data format is [non-standard](https://uk-biobank.gitbook.io/data-access-guide/the-main-dataset/the-structure-of-a-main-dataset) and because some fields contain information about multiple phenotypes. Furthermore, multiple fields can also be combined to define custom phenotypes. The `UKB_CONFIG` file lets you provide a mapping between the UK Biobank Data-Fields and actual phenotypes. By default, it defines 110 non-binary and 660 binary traits as previously defined by the [geneATLAS](http://geneatlas.roslin.ed.ac.uk/). In cases, you may be interested in only a few traits or would like to define more. The goal of the following sections is to help you do that. We first describe some of the main UK-Biobank fields, then explain the rules according to which they are converted into individual phenotypes and end with some examples.

#### Data Fields Description

Phenotypes in the UK-Biobank are organised in Data-Fields, that may contain information about multiple phenotypes or can be combined, here is a non-exhaustive lit of them.

*Diagnoses made during hospital inpatient admission* (467 traits, Data-Fields [41204](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41204) and [41202](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=41202)). This category contains information relating to main and secondary diagnoses made during hospital inpatient admissions. Each trait in this category corresponds to a node or set of nodes in the International Classification of Diseases 10-th Revision (ICD10) ontology. For example, "K41 Femoral hernia" is defined by any of the following diagnoses ("K410", "K411", "K412", "K413", "K414", "K419").

*Self reported conditions* (161 traits, Data-Field [20002](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=20002)). This category contains data obtained through a verbal interview by a trained nurse on past and current medical conditions, including type of cancer and other illnesses, the number of medical conditions, and date of diagnosis. This field can contain inaccuracies due to recall bias or intentional misreporting.

*Self reported conditions* (29 traits, Data-Field [40006](https://biobank.ndph.ox.ac.uk/ukb/field.cgi?id=40006)) This category contains coded data on cancer incidence, obtained through linkage to national cancer registries. Because data is continually accruing, the number of cases may vary depending on the time the data is sent out to researchers.

*Further Miscellaneous fields* (113 traits) These traits do not correspond to diseases but to various lifestyle behaviours or biological measurements. For instance Data-Field 30180 corresponds to Lymphocyte percentage in a blood assay while Data-Field 1389 corresponds to Pork intake.

#### The Extraction Rules

Phenotypes are defined from fields according to rules. These rules are either defined at the field level, or based on its "value type" [metadata](https://biobank.ndph.ox.ac.uk/ukb/schema.cgi?id=1), which indicates the data type the field contains. Here we summarise these rules.

*Continuous and Integer Fields* (value type = 31, 11 respectively). These contain one or more values per individual corresponding to the multiple visit assessments. The value of the field at the first visit assessment is extracted. For instance, the "Alcohol intake frequency." is simply the value from the "1558-0.0" column.

*Ordinal Fields* (1408, 1727, 1548, 728, 1717, 1389, 1478, 1518, 1558, 1349, 1359, 1369, 1379, 1329, 1339, 1239, 1687, 1697, 1319, 1498). Some variables are encoded as categorical by the UK-Biobank, but an ordinal interpretation seems more appropriate. For instance, "Lamb intake" ranges from 0 (Never) to 5 (Once or more daily). Similarly to continuous fields, the value of the field at the first visit assessment is extracted. Negative values (Do not know, Prefer not to answer) are treated as missing.

*Categorical Fields* (value type = 21). These variables cannot be turned into continuous variables and are usually represented via a coding. For instance, the field 21000, described as Ethnic background, uses coding 1001 for British individuals, 4001 for Caribbeans 5 for Chinese etc... These can currently be processed in two different ways by TarGene. Either the coding of each individual is extracted, or, if multiple codings are specified, an indicator of these codings is produced (1 or 0). For instance, "White" could be defined as any of the following codings {1001, 1002, 1003}, that is British, Irish or Any other white background. Similarly to ordinal variables, the value at the first visit assessment is extracted and negative values are considered missing.

*Categorical Arrayed Fields* ({40006, 20002, 41202, 41204}). Some fields in the UK-Biobank comprise a list of codings for each individual at each visit. Each coding represents, for example, a disease or condition that was diagnosed or self reported. A phenotype is then defined by a set of codings. As an example, we define "oesophageal disorder" as the set of the following codings: {1134, 1139, 1140, 1141, 1138, 1474} from Data-Field 20002. An individual annotated with any of these codings, at any assessment visit, is considered a case for "oesophageal disorder". Moreover, some fields share the same codings, this is the case for 41202 and 41204. In that situation it may be useful to aggregate these sources of information. For instance, "G20 Parkinson's disease" is defined by the single element set {G20}, in any of the 41202 and 41204 Data-Fields.

#### The `UKB_CONFIG` Configuration File

We now explain how these rules are encoded within the `UKB_CONFIG` YAML file. For reference, the default [configuration file](https://github.com/TARGENE/targene-pipeline/blob/main/assets/ukbconfig.yaml) corresponds to the [geneATLAS](http://geneatlas.roslin.ed.ac.uk/) study and can be used as a template.

The structure of the `UKB_CONFIG` YAML file consists of extraction rules that convert UK-Biobank fields to traits of interest. It contains two sections, a `traits` section and a `subset` section.

- The `traits` section describes the phenotypes that need to be extracted from the main `TRAITS_DATASET`.
- The `subset` section provides a way to filter the data based on specific traits and hence estimate conditional effects. Each element in this section corresponds to an additional condition to be met by an individual. For instance, the default configuration only includes White individuals.

Each section is further divided in a list of `fields` items that match the UK-Biobank [data fields](https://biobank.ctsu.ox.ac.uk/crystal/help.cgi?cd=data_field) and a list of `phenotypes` that can be extracted from these fields. Each phenotype in the `phenotypes` subsection of the `UKB_CONFIG` YAML file is identified by a `name` and an optional list of `codings` defining it. Altogether, a phenotype element defines an extraction rule from a list of fields and codings.

It is probably more helpful to look at an example:

```yaml
traits:
    - fields:
        - "41202"
        - "41204"
      phenotypes:
        - name: "Chronic lower respiratory diseases"
          codings:
            - "J40"
            - "J41"
            - "J410"
            - "J411"
            - "J418"
            - "J42"
            - "J43"
        - name: "Other extrapyramidal and movement disorders"
          codings:
            - "G254"
            - "G255"
            - "G256"
            - "G258"
    - fields:
        - "40006"
      phenotypes:
        - name: "Malignant melanoma of skin"
          codings:
            - "C430"
            - "C431"
            - "C432"
            - "C433"
            - "C434"
            - "C435"
            - "C436"
            - "C437"
            - "C438"
            - "C439"
    - fields:
        - "1558"
      phenotypes:
        - name: "Alcohol intake frequency."
    - fields:
        - "20117"
      phenotypes:
        - name: "Alcohol drinker status"

subset: 
    - fields: 22001
      phenotypes:
        - name: Female
          codings: 0
    - fields: 21000
      phenotypes:
        - name: White
          codings: [1001]
```

This configuration file restricts the analysis to White females and extract 5 traits.

The first trait, "Chronic lower respiratory diseases" is defined by **any** of the (J40, J41, J410, J411, J418, J42, J43) codings within **any** of the fields (41202, 41204).

The phenotype "Malignant melanoma of skin" is defined by any of (C430, C431, C432, C433, C434, C435, C436, C437, C438, C439) within the 40006 field.

The "Alcohol drinker status" is defined from the first visit assessment questionnaire status.

!!! tip "Constructing Additional columns"
    If you are providing a decrypted `TRAITS_DATASET`, you can build and add extra columns to it (for instance to define new phenotypes). Any column in the `TRAITS_DATASET` which does not correspond to a UK-Biobank field will be considered as such.

#### Additional Files

Additional optional UK-Biobank files for preprocessing and filtering are:

- `QC_FILE`: A path to the UK-Biobank SNP quality control [`ukb_snp_qc.txt`](https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=1955) file.
- `UKB_WITHDRAWAL_LIST`: A path to the withdrawal sample list to exclude removed participants from the study.

### All of Us cohort

TarGene can now be run on the All of Us (AoU) Cohort through the AoU Researcher Workbench. Your Traits Dataset must be built on the AoU Researcher Workbench through their interactive Dataset Builder tool (See [`Dataset Builder`](https://support.researchallofus.org/hc/en-us/articles/4556645124244-Building-a-Dataset-with-the-Dataset-Builder) for more information on this). This must be done on a Workspace created to run your analysis-of-interest (See [`Creating a Workspace`](https://support.researchallofus.org/hc/en-us/articles/30143658322836-Creating-a-Workspace)). The dataset created using these tools must be tailored to pull the traits, confounders and/or covariates-of-interest from the All of Us cohort, subsetting for samples for which genetic data is available. As TarGene requires both PLINK BED and BGEN files, these must be available for the participants for which trait data is being pulled. 

Through Dataset Builder tool on the AoU Researcher Workbench, you must create a matrix that matches the input required by the `COHORT = "CUSTOM"` mode. This requires your partipicant IDs (named, by default, as `person_id` in AoU) to be within a column named `SAMPLE_ID`, and any covariates or confounders to be included in subsequent columns. Each covariate or confounder in the subsequent columns must have a one-to-one mapping with each `SAMPLE_ID`. Please note that the AoU cohort includes data compiled across Electronic Health Care records, and therefore the same patient may have multiple entries for a given measurement (for example, BMI measured across multiple GP appointments). This must be dealt with accordingly when configuring your Traits Dataset. For example, you may choose to pick the most recent measurement for a given participant. The AoU Researcher Workbench provides interactive jupyter notebooks where python can be leveraged to perform these kinds of operations. 

The AoU cohort provides some smaller callsets derived from WGS data, for which both PLINK BED and BGEN files are available (termed `srWGS callsets`). We recommend using these with TarGene if you decide to run analyses using the AoU cohort. For more information about these callsets, see [`Short Read WGS Callsets`](https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK). 

Please be aware that the callset for ACAF-thresholded genetic data is extremely large and may result in longer runtimes. This is because TarGene is run in parallel using the Google Lifesciences API on the AoU Researcher Workbench. Since this is a cloud-based platform, each task is run individually on a runner Virtual Machine (VM), so all data relevant to a task must be copied onto that VM before the task can be executed. This results in a large amount of copying overhead, as symbolic links cannot be leveraged between tasks. 
