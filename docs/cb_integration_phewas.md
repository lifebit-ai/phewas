# CB integration with GWAS pipeline

CB output needs to be prepared and fed into pheWAS pipeline so users have an end-to-end experience.

## 1. Requirements

According to https://lifebit.gitbook.io/cloudos/pipelines-documentations-and-examples-1/nextflow-pipelines/phewas :

### 1.1 - Inputs

- **--vcf_file:** CSV file containing path to VCF files and the following columns:

    1. The `sampleid` column contains a unique sample ID for each individual. These IDs are the same as the IDs used in the file for the `pheno_file` parameter
    2. The `vcf` column contains the path to an input VCFs. Each VCF file contains one individual.
    3. The `PHENO` column contains a phenotype of interest such as SMOKER. This can be used for the `pheno` parameter to select to top SNPs of interest if the snps parameter is not provided. The phenotype values are `1` = control, `2` = case and `-0`/`0`/non-numeric = missing data if case/control
    4. The `SEX` column contains the values `1` = male,  `2` = female, `0` = unknown
    5. Other `PHENO` columns can also be used such as `AGE`

```
sampleid,vcf,SMOKER,SEX,AGE
1,s3://test-vcfs/mini_cohort/1_1.vcf.gz,1,2,39
2,s3://test-vcfs/mini_cohort/2_2.vcf.gz,1,1,52
3,s3://test-vcfs/mini_cohort/3_3.vcf.gz,2,1,42
4,s3://test-vcfs/mini_cohort/4_4.vcf.gz,2,2,44
5,s3://test-vcfs/mini_cohort/5_5.vcf.gz,2,2,23
```

OR

- **.bed .bim .fam files**


- **PhenoCol** Name of the column in the .fam/VCF_file (i.e SMOKER)
- **PhenoFile:** Contains DOID to be mapped to ICD9 or ICD10 codes and counts
```
id,doid,count
1170,3055,4
1170,12750,4
1170,12385,4
1170,1498,4
1170,12385,9
1170,3055,10
1170,12385,5
1170,0060859,3
1170,12385,3
1170,14239,3
```
### 1.2 - Outputs

### 1.3 - Scripts


# Future work

- [ ] **(1) Automatic inference of sex from genotypic data**
   - [ ] Check how to do it: 
      - seXY: https://academic.oup.com/bioinformatics/article/33/4/561/2666346
      - plink: https://www.cog-genomics.org/plink/1.9/basic_stats#check_sex 
   - [ ] Create trigget to run this option when no sex in phenoFile or user ask for it (exposing parameter).
   - [ ] Make a Nextflow process which runs sex inference and adds it to the phenoFile

- [ ] **(2) Branch GWAS processes so that multiple reports can be generated and explored:**
   - [ ] Modify reporting process so it allows to combine multiple reports in one report
   - [ ] Modify reporting process so it allows to switch between reports

- [ ] **(3) GWAS post-analysis**
   - [ ] Add MAGMA & friends to run pathway-analysis
   - [ ] Add IntAssoPlot https://www.frontiersin.org/articles/10.3389/fgene.2020.00260/full
   - [ ] PRS within pipeline