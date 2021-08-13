# General Usage

In order to use this pipeline, you can run the following example:

```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10"

```
`input_phenofile` is only required if `--snps` parameter is **NOT** used. The phenotype is used for association testing to derive at a list of SNPs (most significant hits) to use in the pheWAS step. 

If you are specifying a list of SNPs via `--snps` parameter, `--input_samplefile` can be used instead of `--input_phenofile`.

## Parameters

### **ESSENTIAL**

#### Genomic data

- **--plink_input** : Path/URL to plink bim/bed/fam files e.g. /path/plink{.bed,.bim,.fam}

OR individual plink files can be supplied:

- **--bim** : Path/URL to .bim file.
- **--bed** : Path/URL to .bed file.
- **--fam** : Path/URL to .fam file.

OR:

- **--agg_vcf_file_list** : Path/URL to a .csv file containing chr chunk information, path to aggregated VCFs, VCFs index. The required column headers are `chr`, `vcf` and `index`, e.g.
```
chr,vcf,index
1,s3://test-vcfs/mini_cohort/chr1.vcf.gz,1,s3://test-vcfs/mini_cohort/chr1.vcf.gz.csi
2,s3://test-vcfs/mini_cohort/chr2.vcf.gz,s3://test-vcfs/mini_cohort/chr2.vcf.gz.csi
3,s3://test-vcfs/mini_cohort/chr3.vcf.gz,s3://test-vcfs/mini_cohort/chr3.vcf.gz.csi
4,s3://test-vcfs/mini_cohort/chr4.vcf.gz,s3://test-vcfs/mini_cohort/chr4.vcf.gz.csi
5,s3://test-vcfs/mini_cohort/chr5.vcf.gz,s3://test-vcfs/mini_cohort/chr5.vcf.gz.csi
```


OR:

- **--individual_vcf_file_list** : Path/URL to .csv file containing individual, path to individual VCFs. The required column headers are `sampleid` and `vcf`, e.g.
```
sampleid,vcf
1_1,s3://test-vcfs/mini_cohort/1_1.vcf.gz
2_2,s3://test-vcfs/mini_cohort/2_2.vcf.gz
3_3,s3://test-vcfs/mini_cohort/3_3.vcf.gz
4_4,s3://test-vcfs/mini_cohort/4_4.vcf.gz
5_5,s3://test-vcfs/mini_cohort/5_5.vcf.gz
6_6,s3://test-vcfs/mini_cohort/6_6.vcf.gz
7_7,s3://test-vcfs/mini_cohort/7_7.vcf.gz
8_8,s3://test-vcfs/mini_cohort/8_8.vcf.gz
```


#### Phenotypic data

- **--input_phenofile** : Path/URL to file containing phenotypic information about input cohort. Only required if `--snps` parameter is used, otherwise `--input_samplefile` can be used instead.
- **--input_id_code_count** : Path/URL to file that contains patient id, pheno code, counts from phenotypic information, e.g. 
```
id,code,count
1_1,A01.0,3
2_2,A15,3
2_2,A06.1,3
3_3,A01.0,3
3_3,A06.1,3
```
The required columns are `id`, `code` and `count` where `id` contains sample IDs corresponding to genotype data, `code` contains supplied disease classification groups (corresponding to the choice of `--pheno_codes`) and `count` refers to counts of each code per individual.

- **--pheno_codes** : String containing phenotypic codes to be used. Select from "icd9", "doid", "hpo" or "icd10". Default = "doid".

## **Optional**

- **--concat_vcfs** Combines input VCF files by concatenation (bcftools concat) - recommended for per-chromosome inputs. Default behaviour is `bcftools merge`. `--concat_vcfs` should be specified as a flag, with no input.

- **--snps** File containing relevant SNPs to be tested in a `phewas` step. The file should contain no header and one SNP per row. The variant IDs should correspond to those found in VCF or .bim files.
```
rs1002365
rs1023115
rs1041676
rs1048418
rs10489265
rs10493582
rs10493644
rs10494452
rs10494464
```

- **--snp_p_val_threshold** : Threshold to select SNPs significance from plink association analysis if `snps` are not provided. Default = 0.05.

- **--outdir** : Output directory for results.
- **--output_tag** : Prefix to identify output files.
- **--post_analysis** : String containing type of posterior analysis to be executed. Only option is `"coloc"`
- **--gwas_input** : Path to file containing GWAS summary statistics.
- **--gwas_trait_type** : String containing type of trait in GWAS. Accepts `"binary"` or `"quantitative"`.

## **Binary**

Only binary GWAS type is currently implemented.

## **Quantitative**

Option not implemented or tested yet.

## **Colocalization**

**Binary GWAS**

```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis 'coloc' \
                     --gwas_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/gwas_summary_bin.csv" \
                     --gwas_trait_type 'binary'

```

**Quantitative GWAS**

```bash
nextflow run main.nf --input_phenofile "s3://lifebit-featured-datasets/pipelines/phewas/cb_binary_pheno.phe" \
                     --input_id_code_count "s3://lifebit-featured-datasets/pipelines/phewas/phewas_id_code_count.csv" \
                     --plink_input "s3://lifebit-featured-datasets/pipelines/biobank-gwas/testdata/sampleA.{bed,bim,fam}" \
                     --pheno_codes "icd10" \
                     --post_analysis 'coloc' \
                     --gwas_input 's3://lifebit-featured-datasets/pipelines/biobank-gwas/gwas_summary_qt.csv' \
                     --gwas_trait_type 'quantitative'

```

### 1.2 - Outputs
`phewas` folder containing:
- `phewas_man.png` PheWAS Manhattan plot

`plink` folder containing:
- `r_genotypes.raw` contains information from the FAM file with the number of alleles for each given SNP of interest for all individuals (can be loaded into R)
- `r_genotypes.log` log file detailing the execution of the PLINK command

`phewas` folder containing results (as many files as chromosomes or one unique file) from pheWAS:
- `*_phewas_results.csv` summary statistics.


`MultiQC` directory containing:
- `.report.json` file used to render the visualisations within CloudOS
