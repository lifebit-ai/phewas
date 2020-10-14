# Usage

In order to use this pipeline, you can run the following example:
**Binary**
```bash
nextflow run main.nf --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/76420a552c7f3bae7619fc2d56605ad06165ea84/cohort_data_phenos_phewas.csv" \
                     --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
                     --continuous_var_transformation "mean" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --phenoCol "Specimen type" \
                     --case_group "NOSE" \
                     --mode "case_vs_control_contrast" \
                     --pheno_codes "icd10"
```

# Parameters

## **ESSENTIAL**

- **--plink_input** : path/url to CSV file containing chr chunk information, path to aggregated VCFs, VCFs index.
- **--phenoCol** : String with the name of the phenotypic column to be used as trait. Note for CB users, it must match the name of the column, i.e. 'Specimen type'.
- **--phenofile** : path/url to file that contains phenotypic information about cohort to be analysed.
- **--metadata** : path/url to file that contains metadata from phenotypic information.
- **--pheno_codes** : String containing phenotypic codes to be used. Select from "icd9", "doid", "hpo or "icd10".

## **Optional**

- **continuous_var_transformation** : Transforms continuous variables using 'log', 'log10', 'log2', 'zscores' or 'None'.
- **continuous_var_aggregation** : Defines how to aggregate different measurements. Choose between 'max', 'min', 'median' or 'mean'.
- **snps** : File containing list of SNPs to be tested.
- **snp_threshold** : Threshold to select SNPs significance from plink association analysis if `snps` are not provided. Defaults to 0.05
- **outdir** : Output directory for results.
- **output_tag** : Prefix to identify output files.

## **Binary**
- **--mode** : String containing the design matrix configuration to be used. Allows the user to select between the following scenarios:

|| Value | Description | Needs | Added value |
|--|--|--|--|--|
| Scenario 1 | case_vs_control_contrast | User wants to run on a particular case but wants to use all the rest of cases as controls. | Subset the particular case group and select all the remaining individuals as control | Find significant associations exclusive to the case group you are interested |
| Scenario 2 | case_vs_group_contrasts | User wants to run on a particular case but wants to compare to each of the other cases as controls independently | Subset case vs each group as control | Find associations that are different to an specific group |
| Scenario 3 | all_vs_all | User doesn't have a particular group in mind and wants to run an exploration on the phenotype | All vs All approach | Allows for exploration or assumptions free analysis |

- **case_group** : String containing name of the case group selected for contrasts.

## **Quantitative**

Option not implemented and tested yet.
