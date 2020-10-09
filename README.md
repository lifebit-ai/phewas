# phewas

Example commands to run the pipeline:
```bash
nextflow run main.nf --vcf_file testdata/vcf_samples.csv --pheno SMOKER --pheno_file testdata/all.doid.count.csv
```
Or with CB outputs:
```bash
nextflow run main.nf --phenofile "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/76420a552c7f3bae7619fc2d56605ad06165ea84/cohort_data_phenos_phewas.csv" \
                     --metadata "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/d602bec4b31d5d75f74f1dbb408bd392db57bdb6/metadata.csv" \
                     --continuous_var_transformation "mean" \
                     --plink_input "testdata/plink/sampleA_*.{bed,bim,fam}" \
                     --phenoCol "Specimen type" \
                     --case_group "NOSE" \
                     --mode "case_vs_control_contrast" \
                     --pheno_codes "icd10"
```