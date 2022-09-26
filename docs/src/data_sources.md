# Setting a data source

Currently, only the UK-Biobank is supported, stay tuned for further data sources!

## UK-Biobank

The UK-Biobank is composed of both genetic data (.bed and .bgen files) and trait data. While there is a plan to use only bgen and trait data, the pipeline currently uses all sources of information. Below we explain how to specify those arguments to TarGene. For more information on the structure of the UK-Biobank data, please refer to their [User Guide](https://biobank.ndph.ox.ac.uk/crystal/exinfo.cgi?src=accessing_data_guide).

### Main Dataset

This dataset is often called a "main dataset" and consists of an encrypted file containing individuals trait information accessible via your project. The first option is thus to provide this dataset using the `ENCRYPTED_DATASET` parameter. Since the data is encrypted, the pipeline will also need the encoding file that you can provide with `ENCODING_FILE`.

A "main dataset" is typically very large and only a few traits will be of interest to a given study. To extract those relevant traits from the dataset, a `TRAITS_CONFIG` YAML file must be provided. The structure of a file is composed of extraction rules that convert UK-Biobank fields to traits of interest. An example is presented below:

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

A file consists in two sections, a `traits` section and a `subset` section. Each section further consists in a list of `fields` items. Each `fields` item is itself decomposed in multiple `phenotypes` items that are identified by a `name` and an optional list of `codings`. Altogether, a phenotype element defines an extraction rule. The previous configuration file would thus restrict the analysis to white females and extract 5 traits. Since, writing by hand such a file for large scale study can quickly become tenuous, we provide a configuration file corresponding to the [GeneAtlas](http://geneatlas.roslin.ed.ac.uk/) study [here]().

We now describe the currently available extraction rules based on variable types. Not that a variable type is determined by the [metadata](https://biobank.ndph.ox.ac.uk/ukb/schema.cgi?id=1) provided by the UK-Biobank and not a programmatic rule.

- Continuous/Count variables: Those consist in a single `field` and the extraction rule will take the first visit assessment value for that field. For instance, the "Alcohol intake frequency." is simply the value from the "1558-0.0" column.
- Ordinal data: Some variables are encoded as categorical by the UK-Biobank but an ordinal interpretation seems more appropriate. A list of such variables has been identified and hard-coded. They will be treated like continuous variables.
- Categorical variables: For those variables again, the first visit assessment value for the specified field is used. An additional `codings` field can be provided. In that case the variable is turned into a indicator variable indicating whether an individual is included in the criterion given by the `codings` list. In the previous `subset` section, "White" is such a transformation into a binary variable. Note that we could add more elements to the list to subset on a wider population.
- Categorical Arrayed variables: Some fields in the UK-Biobank consist in list of traits for individuals. This is the case for at least the fields: 41202, 41204, 20002 and 40006. For those fields it is essential to define a `codings` section that describes a trait. For instance "Malignant melanoma of skin" corresponds to the declaration of any of the `codings` for an individual. Moreover, some fields share the same encoding, this is the case for  41202 and 41204. In that situation it may be useful to aggregate those sources of information. For that purpose multiple fields can be provided to the `fields` section. In that scenario the declaration of any of the coding in a `codings` section for any of the `fields` will define the trait.

To come back to the pipeline specification, if a "typical" main dataset has already been decrypted outside of the pipeline, one may use the `DECRYPTED_DATASET` as an input to the pipeline instead of the `ENCRYPTED_DATASET` and `ENCODING_FILE` arguments. However, the `TRAITS_CONFIG` file described above is still mandatory.

### Genetic Data

We are currently using both .bgen and .bed files. Those are respectively provided with the `UKBB_BGEN_FILES` and `UKBB_BED_FILES` parameters. Since the UK-Biobank genotypic data is split in chromosomes, it should be of the form `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}` and `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}` respectively.

### Additional Files

Additional UK-Biobank required files for preprocessing and filtering are:

- `QC_FILE`: A path to the UK-Biobank SNP quaility control `ukb_snp_qc.txt` file.
- `WITHDRAWAL_LIST`: A path to the withdrawal sample list to exclude removed participants from the study.
