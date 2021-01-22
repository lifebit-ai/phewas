# General Usage

In order to use this pipeline, you can run the following example:

**Binary**

```bash
nextflow run main.nf --input_phenofile "s3://marcos-lifebit/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://marcos-lifebit/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10"

```

# Parameters

## **ESSENTIAL**

- **--plink_input** : path/url to CSV file containing chr chunk information, path to aggregated VCFs, VCFs index.
- **--input_phenofile** : path/url to file that contains phenotypic information about cohort to be analysed.
- **--input_id_code_count** : path/url to file that contains patient id, pheno code, counts from phenotypic information.
- **--pheno_codes** : String containing phenotypic codes to be used. Select from "icd9", "doid", "hpo or "icd10".

## **Optional**

- **--snps** : File containing list of SNPs to be tested.
- **--snp_p_val_threshold** : Threshold to select SNPs significance from plink association analysis if `snps` are not provided. Defaults to 0.05

- **--outdir** : Output directory for results.
- **--output_tag** : Prefix to identify output files.
- **--post_analysis** : String containing type of posterior analysis to be executed. Only option is `"coloc"`
- **--gwas_input** : Path to file containing GWAS summary statistics.
- **--gwas_trait_type** : String containing type of trait in GWAS. Accepts `"binary"` or `"quantitative"`.

## **Binary**

Only accepts this kind of traits.

## **Quantitative**

Option not implemented or tested yet.

## **Colocalization**

**Binary GWAS**

```bash
nextflow run main.nf --input_phenofile "s3://marcos-lifebit/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://marcos-lifebit/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis 'coloc' \
                     --gwas_input "https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/74e0e3b0f1a9c5f95804053b375258da3bfe64cc/gwas_summary_bin.csv" \
                     --gwas_trait_type 'binary'

```

**Quantitative GWAS**

```bash
nextflow run main.nf --input_phenofile "s3://marcos-lifebit/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://marcos-lifebit/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/projects/gel/gel-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis 'coloc' \
                     --gwas_input 'https://gist.githubusercontent.com/mcamarad/e98cdd5e69413fb6189ed70405c43ef4/raw/74e0e3b0f1a9c5f95804053b375258da3bfe64cc/gwas_summary_qt.csv' \
                     --gwas_trait_type 'quantitative'

```
