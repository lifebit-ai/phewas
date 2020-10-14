# phewas

Example commands to run the pipeline:
```bash
nextflow run main.nf --vcf_file testdata/vcf_samples.csv --pheno SMOKER --pheno_file testdata/all.doid.count.csv
```
Or with CB outputs:
```bash
nextflow run main.nf --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/76420a552c7f3bae7619fc2d56605ad06165ea84/cohort_data_phenos_phewas.csv" \
                     --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
                     --continuous_var_aggregation "mean" \
                     --continuous_var_transformation "log10" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_col "Specimen type" \
                     --case_group "NOSE" \
                     --design_mode "case_vs_control_contrast" \
                     --pheno_codes "icd10"
```
With colocalization analysis with binary GWAS:
```bash
nextflow run main.nf --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/76420a552c7f3bae7619fc2d56605ad06165ea84/cohort_data_phenos_phewas.csv" \
                     --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
                     --continuous_var_aggregation "mean" \
                     --continuous_var_transformation "log10" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_col "Specimen type" \
                     --case_group "NOSE" \
                     --design_mode "case_vs_control_contrast" \
                     --pheno_codes "icd10" \
                     --post_analysis "coloc" \
                     --gwas_input "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/74e0e3b0f1a9c5f95804053b375258da3bfe64cc/gwas_summary_bin.csv" \
                     --gwas_trait_type "binary"

```
With colocalization analysis with quantitative GWAS:
```bash
nextflow run main.nf --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/76420a552c7f3bae7619fc2d56605ad06165ea84/cohort_data_phenos_phewas.csv" \
                     --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
                     --continuous_var_aggregation "mean" \
                     --continuous_var_transformation "log10" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_col "Specimen type" \
                     --case_group "NOSE" \
                     --design_mode "case_vs_control_contrast" \
                     --pheno_codes "icd10" \
                     --post_analysis "coloc" \
                     --gwas_input "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/74e0e3b0f1a9c5f95804053b375258da3bfe64cc/gwas_summary_qt.csv" \
                     --gwas_trait_type "quantitative"

```