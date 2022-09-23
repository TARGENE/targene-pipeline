# Setting a data source

Currently, only the UK-Biobank is supported, stay tuned for further data sources!

## UK-Biobank

The UK-Biobank is composed of both genetic data (.bed and .bgen files) and trait data.

- The trait data: The first option is to provide the path to the encrypted dataset `ENCRYPTED_DATASET` which must be decoded via the ukbconv software and the encoding file `ENCODING_FILE`. The second option is to decrypt the dataset only once outside of the pipeline and then use the `DECRYPTED_DATASET` as an input to the pipeline. Finally, since one is usually not interested in all of the traits, and those traits can play a different role in the causal model, the `TRAITS_CONFIG` file described [here](https://github.com/TARGENE/UKBMain.jl) has to be provided.

- The genetic data: We are currently using both .bgen and .bed files. Those are respectively provided with the `UKBB_BGEN_FILES` and `UKBB_BED_FILES` parameters. Since the UK-Biobank genotypic data is split in chromosomes, it should be of the form `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bgen,sample,bgen.bgi}` and `PREFIX_TO_CHROMOSOMES{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}.{bed,bim,fam}` respectively.

Additional UK-Biobank required files for preprocessing and filtering are:

- `QC_FILE`: A path to the UK-Biobank SNP quaility control `ukb_snp_qc.txt` file.
- `WITHDRAWAL_LIST`: A path to the withdrawal sample list to exclude removed participants from the study.