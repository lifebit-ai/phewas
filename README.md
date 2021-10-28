# phewas

Example commands to run the pipeline:
```bash
nextflow run main.nf --agg_vcf_file_list testdata/vcf_samples.csv \
                     --input_phenofile testdata/cb_binary_pheno.phe \
                     --input_id_code_count testdata/phewas_id_code_count.csv \
                     --pheno_codes "icd10"
```
Alternatively, it is possible to specify plink input files:
```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10"
```
A covariate file can be used for pheWAS analysis (optional, but strongly recommended). Please note that the format for the covariate file should be a comma-delimited file with sample IDs in first column (`id` as first column name). The regression model used in pheWAS will use all supplied covariates by default :
```
id,sex,age...
1,1,25
```

```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --covariate_file "s3://eu-west-2-example-data/pipelines/phewas/testdata/cov.csv"
                     --pheno_codes "icd10"
```


`--firth_regression` parameter can supplied to run the association tests with Firth logistic regression - this is highly recommended for phenotypes with significantly imbalanced case/control ratios.


With colocalization analysis with binary GWAS:
```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis "coloc" \
                     --gwas_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/gwas_summary_bin.csv" \
                     --gwas_trait_type "binary"

```
With colocalization analysis with quantitative GWAS:
```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis "coloc" \
                     --gwas_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/gwas_summary_qt.csv" \
                     --gwas_trait_type "quantitative"

```
