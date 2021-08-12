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