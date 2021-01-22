# phewas

Example commands to run the pipeline:
```bash
nextflow run main.nf --vcf_file testdata/vcf_samples.csv --pheno SMOKER --pheno_file testdata/all.doid.count.csv
```
Or with CB outputs:
```bash
nextflow run main.nf --input_phenofile "s3://marcos-lifebit/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://marcos-lifebit/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10"
```
With colocalization analysis with binary GWAS:
```bash
nextflow run main.nf --input_phenofile "s3://marcos-lifebit/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://marcos-lifebit/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis "coloc" \
                     --gwas_input "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/74e0e3b0f1a9c5f95804053b375258da3bfe64cc/gwas_summary_bin.csv" \
                     --gwas_trait_type "binary"

```
With colocalization analysis with quantitative GWAS:
```bash
nextflow run main.nf --input_phenofile "s3://marcos-lifebit/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://marcos-lifebit/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis "coloc" \
                     --gwas_input "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/74e0e3b0f1a9c5f95804053b375258da3bfe64cc/gwas_summary_qt.csv" \
                     --gwas_trait_type "quantitative"

```